#include <ctime>
#include <cstdio>           // FILE, popen, feof, fgets
#include <sstream>
#include <fstream>
#include "util.h"


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

string git::commit()
{
    return util::get_command_output("git rev-parse HEAD");
}

bool git::has_uncommitted_changes()
{
    string changes = util::get_command_output(
            "git status --untracked-files=no --porcelain");
    return !changes.empty();
}

string git::commit_info()
{
    try {
        string info = git::commit();
        if (git::has_uncommitted_changes())
            info += " (uncommitted changes)";
        return info;
    }
    catch (std::runtime_error& e)
    {
        return "(failed to get git commit)";
    }
}


void terminal::cursor_up(std::ostream& out, int num_lines)
{
    out << "\033[" << num_lines << "A" << std::flush;
}


void terminal::clear_current_line(std::ostream& out)
{
    out << "\r\033[K" << std::flush;
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
