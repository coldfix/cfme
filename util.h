#ifndef __UTIL_H__INCLUDED__
#define __UTIL_H__INCLUDED__

# include <ctime>       // time_t
# include <iostream>
# include <sstream>
# include <string>
# include <utility>     // move
# include <vector>

# include <boost/timer/timer.hpp>

# include <termios.h>


namespace terminal
{
    void cursor_up(std::ostream&, int num_lines=1);
    void clear_current_line(std::ostream& out);

    class Input
    {
        termios ttystate_restore, ttystate_replace;
        int fd;
    public:
        Input();
        ~Input();

        bool avail();
        int get();
    };
}


namespace util
{
    class AutogenNotice
    {
        boost::timer::cpu_timer timer;
        std::vector<std::string> argv;
        std::time_t start_time;

    public:
        AutogenNotice(int argc, char** argv);

        std::string str() const;
    };

    std::string get_command_output(const std::string& command);

    std::string join(const std::vector<std::string>&, const std::string& sep);

    std::vector<std::string> read_file(std::istream&);
    std::vector<std::string> read_file(const std::string& filename);


    inline void print_all(std::ostream& out)
    {
        out << std::flush;
    }

    template <class H, class ...T>
    void print_all(std::ostream& out, H head, T... t)
    {
        out << head;
        print_all(out, std::move(t)...);
    }

    template <class ...T>
    std::string sprint_all(T... t)
    {
        std::ostringstream out;
        print_all(out, std::move(t)...);
        return out.str();
    }


    template <class T>
    void extend(std::vector<T>& a, const std::vector<T>& b)
    {
        a.reserve(a.size() + b.size());
        a.insert(a.end(), b.begin(), b.end());
    }
}


namespace git
{
    extern const char* commit;
    extern bool has_uncommitted_changes;
    std::string commit_info();
}

#endif // include guard
