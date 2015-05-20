#include <ctime>
#include <cstdio>           // FILE, popen, feof, fgets
#include <sstream>
#include <fstream>
#include "util.h"

#include <stdlib.h>         // these are for terminal Input
#include <termios.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>       // open
#include <fcntl.h>          // O_RDONLY


using std::string;


string util::get_command_output(const string& command)
{
    std::FILE* pipe = popen(command.c_str(), "r");
    if (!pipe)
        throw std::runtime_error("popen failed");

    char buffer[256];
    string result = "";
    while (!std::feof(pipe)) {
    	if (std::fgets(buffer, sizeof(buffer), pipe) != NULL)
    		result += buffer;
    }
    pclose(pipe);
    if (!result.empty() && result.back() == '\n')
        result.pop_back();
    return result;
}


string git::commit_info()
{
    string info = git::commit;
    if (git::has_uncommitted_changes)
        info += " (uncommitted changes)";
    return info;
}


void terminal::cursor_up(std::ostream& out, int num_lines)
{
    out << "\033[" << num_lines << "A" << std::flush;
}


void terminal::clear_current_line(std::ostream& out)
{
    out << "\r\033[K" << std::flush;
}


terminal::Input::Input()
{
    fd = open("/dev/tty", O_RDONLY);
    tcgetattr(fd, &ttystate_restore);
    ttystate_replace = ttystate_restore;
    ttystate_replace.c_lflag &= ~(ICANON | ECHO);
    ttystate_replace.c_cc[VMIN] = 1;        // minimum of number input read.
    tcsetattr(fd, TCSANOW, &ttystate_replace);
}

terminal::Input::~Input()
{
    tcsetattr(fd, TCSANOW, &ttystate_restore);
    close(fd);
}

int terminal::Input::get()
{
    char c;
    read(fd, &c, 1);
    return c;
}

bool terminal::Input::avail()
{
    // from http://cc.byexamples.com/2007/04/08/non-blocking-user-input-in-loop-without-ncurses/
    struct timeval tv;
    fd_set fds;
    tv.tv_sec = 0;
    tv.tv_usec = 0;
    FD_ZERO(&fds);
    FD_SET(fd, &fds);
    select(fd+1, &fds, NULL, NULL, &tv);
    return FD_ISSET(fd, &fds);
}


string util::join(const std::vector<string>& v, const string& sep)
{
    if (v.empty())
        return std::string();
    std::string result = v.front();
    for (int i = 1; i < v.size(); ++i) {
        result += sep;
        result += v[i];
    }
    return result;
}

std::vector<std::string> util::read_file(std::istream& in)
{
    std::vector<string> lines;
    string line;
    while (std::getline(in, line))
        lines.push_back(line);
    return lines;
}

std::vector<string> util::read_file(const string& filename)
{
    std::ifstream in(filename);
    return read_file(in);
}

util::AutogenNotice::AutogenNotice(int _argc, char** _argv)
    : argv(_argv, _argv + _argc)
    , start_time(std::time(nullptr))
{
}

string util::AutogenNotice::str() const
{
    using std::endl;
    std::ostringstream out;
    out << "# command line: " << util::join(argv, " ") << endl;
    out << "# start date:   " << std::ctime(&start_time);
    out << "# git commit:   " << git::commit_info() << endl;
    out << "# running time: " << timer.format(3) << endl;
    string result = out.str();
    result.pop_back();
    return result;
}
