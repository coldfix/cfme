#ifndef __UTIL_H__INCLUDED__
#define __UTIL_H__INCLUDED__

#include <iostream>
#include <string>
#include <vector>

#include <boost/timer/timer.hpp>
#include <ctime>


namespace terminal
{
    void clear_current_line(std::ostream& out);
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
}


namespace git
{
    std::string commit();
    bool has_uncommitted_changes();
    std::string commit_info();
}

#endif // include guard
