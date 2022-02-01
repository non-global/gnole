///////////////////////////////////////////////////////////////////////////////
// File: CmdLine.cc                                                          //
// Part of the CmdLine library                                               //
//                                                                           //
// Copyright (c) 2007-2032 Gavin Salam with contributions from               //
// Gregory Soyez and Rob Verheyen                                            //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "CmdLine.hh"
#include<string>
#include<sstream>
#include<iostream> // testing
#include<fstream>
#include<vector>
#include<cstddef> // for size_t
#include <sys/utsname.h> // for getting uname
#include <unistd.h> // for getting current path
#include <stdlib.h> // for getting the environment (including username)
#include <cstdio>
using namespace std;

string CmdLine::_default_argfile_option = "-argfile";

// initialise the various structures that we shall
// use to access the command-line options;
//
// If an option appears several times, it is its LAST value
// that will be used in searching for option values (opposite of f90)
CmdLine::CmdLine (const int argc, char** argv, bool enable_help, const string & file_option) : 
    __help_enabled(enable_help), __argfile_option(file_option) {

  __arguments.resize(argc);
  for(int iarg = 0; iarg < argc; iarg++){
    __arguments[iarg] = argv[iarg];
  }
  this->init();
}

/// constructor from a vector of strings, one argument per string
CmdLine::CmdLine (const vector<string> & args, bool enable_help, const string & file_option) : 
    __help_enabled(enable_help), __argfile_option(file_option) {

  __arguments = args;
  this->init();
}

/// Add an overall help string
CmdLine & CmdLine::help(const std::string & help_str) {
  __overall_help_string = help_str;
  __help_enabled = true;
  return *this;
}

//----------------------------------------------------------------------
void CmdLine::init (){
  // record time at start
  time(&__time_at_start);

  // this does not work...
  //__options_help[__argfile_option] = 
  //    OptionHelp_value_with_default<string>(__argfile_option, "filename", 
  //                          "if present, further arguments are read from the filename");

  // check first if a file option is passed
  for(size_t iarg = 0; iarg < __arguments.size(); iarg++) {
    const string & arg = __arguments[iarg];
    if (arg == __argfile_option) {
      // make sure a file is passed too
      bool found_file = true;
      ifstream file_in;
      if (iarg+1 == __arguments.size()) found_file = false;
      else {
        file_in.open(__arguments[iarg+1].c_str());
        found_file = file_in.good();
      }

      // error if no file found
      if (!found_file) {
        ostringstream ostr;
        ostr << "Option "<< __argfile_option
             <<" is passed but no file was found"<<endl;
        throw Error(ostr);
      }

      // remove the file options from the list of arguments
      __arguments.erase(__arguments.begin()+iarg, __arguments.begin()+iarg+2);

      string read_string = "";
      while (file_in >> read_string) {
        // skip the rest of the line if it's a comment
        if (read_string.find("//") != string::npos) {
          getline(file_in, read_string);
        }
        else {
          __arguments.push_back(read_string);
        }
      }

      // start from the beginning of the argument list again again
      iarg = 0;
    }
  }

  // record whole command line so that it can be easily reused
  __command_line = "";
  for(size_t iarg = 0; iarg < __arguments.size(); iarg++){
    const string & arg = __arguments[iarg];
    // if an argument contains special characters, enclose it in
    // single quotes [NB: does not work if it contains a single quote
    // itself: treated below]
    if (arg.find(' ') != string::npos ||
        arg.find('|') != string::npos ||
        arg.find('<') != string::npos || 
        arg.find('>') != string::npos || 
        arg.find('"') != string::npos || 
        arg.find('#') != string::npos) {
      __command_line += "'"+arg+"'";
    } else if (arg.find("'") != string::npos) {
      // handle the case with single quotes in the argument
      // (NB: if there are single and double quotes, we are in trouble...)
      __command_line += '"'+arg+'"';
    } else {
      __command_line += arg;
    }
    __command_line += " ";
  }
  
  // group things into options
  bool next_may_be_val = false;
  string currentopt;
  for(size_t iarg = 1; iarg < __arguments.size(); iarg++){
    // if expecting an option value, then take it (even if
    // it is actually next option...)
    if (next_may_be_val) {__options[currentopt] = iarg;}
    // now see if it might be an option itself
    string arg = __arguments[iarg];
    bool thisisopt = (arg.compare(0,1,"-") == 0);
    if (thisisopt) {
      // set option to a standard undefined value and say that 
      // we expect (possibly) a value on next round
      currentopt = arg;
      __options[currentopt] = -1;
      __options_used[currentopt] = false;
      next_may_be_val = true;}
    else {
      // otherwise throw away the argument for now...
      next_may_be_val = false;
      currentopt = "";
    }
  }
  if (__help_enabled) {
    __help_requested = present("-h").help("prints this help message") || present("--help");
  }
}

// indicates whether an option is present
CmdLine::Result<bool> CmdLine::present(const string & opt) const {
  OptionHelp * opthelp = 0;
  if (__help_enabled && __options_help.find(opt) == __options_help.end()) {
    __options_queried.push_back(opt);
    __options_help[opt] = OptionHelp_present(opt);
    opthelp = & __options_help[opt];
  }
  bool result = (__options.find(opt) != __options.end());
  if (result) __options_used[opt] = true;
  return Result<bool>(result,opthelp);
}

// indicates whether an option is present and has a value associated
bool CmdLine::present_and_set(const string & opt) const {
  bool result = present(opt) && __options[opt] > 0;
  return result;
}


// return the string value corresponding to the specified option
string CmdLine::string_val(const string & opt) const {
  if (!this->present_and_set(opt)) {
    ostringstream ostr;
    ostr << "Option "<<opt
	 <<" is needed but is not present_and_set"<<endl;
    throw Error(ostr);
  }
  string arg = __arguments[__options[opt]];
  // this may itself look like an option -- if that is the case
  // declare the option to have been used
  if (arg.compare(0,1,"-") == 0) {__options_used[arg] = true;}
  return arg;
}

// as above, but if opt is not present_and_set, return default
string CmdLine::string_val(const string & opt, const string & defval) const {
  if (this->present_and_set(opt)) {return string_val(opt);} 
  else {return defval;}
}

// Return the integer value corresponding to the specified option;
// Not too sure what happens if option is present_and_set but does not
// have string value...
int CmdLine::int_val(const string & opt) const {
  int result;
  string optstring = string_val(opt);
  istringstream optstream(optstring);
  optstream >> result;
  if (optstream.fail()) {
    ostringstream ostr;
    ostr << "could not convert option ("<<opt<<") value ("
	 <<optstring<<") to int"<<endl; 
    throw Error(ostr);
  }
  return result;
}

// as above, but if opt is not present_and_set, return default
int CmdLine::int_val(const string & opt, const int & defval) const {
  if (this->present_and_set(opt)) {return int_val(opt);} 
  else {return defval;}
}


// Return the integer value corresponding to the specified option;
// Not too sure what happens if option is present_and_set but does not
// have string value...
double CmdLine::double_val(const string & opt) const {
  double result;
  string optstring = string_val(opt);
  istringstream optstream(optstring);
  optstream >> result;
  if (optstream.fail()) {
    ostringstream ostr;

    ostr << "could not convert option ("<<opt<<") value ("
	 <<optstring<<") to double"<<endl; 
    throw Error(ostr);
  }
  return result;
}

// as above, but if opt is not present_and_set, return default
double CmdLine::double_val(const string & opt, const double & defval) const {
  if (this->present_and_set(opt)) {return double_val(opt);} 
  else {return defval;}
}


// return the full command line including the command itself
string CmdLine::command_line() const {
  return __command_line;
}


// return true if all options have been asked for at some point or other
bool CmdLine::all_options_used() const {
  bool result = true;
  for(map<string,bool>::const_iterator opt = __options_used.begin();
      opt != __options_used.end(); opt++) {
    bool this_one = opt->second;
    if (! this_one) {cerr << "Option "<<opt->first<<" unused/unrecognized"<<endl;}
    result = result && this_one;
  }
  return result;
}

/// return a time stamp corresponding to now
string CmdLine::time_stamp(bool utc) const {
  time_t timenow;
  time(&timenow);
  return _string_time(timenow, utc);
}

/// return a time stamp corresponding to start time
string CmdLine::time_stamp_at_start(bool utc) const {
  return _string_time(__time_at_start, utc);
}

/// return the elapsed time in seconds since the CmdLine object was
/// created
double CmdLine::time_elapsed_since_start() const {
  time_t timenow;
  time(&timenow);
  return std::difftime(timenow, __time_at_start);
}


/// convert the time into a string (local by default -- utc if 
/// utc=true).
string CmdLine::_string_time(const time_t & time, bool utc) const {
  struct tm * timeinfo;
  if (utc) {
    timeinfo = gmtime(&time);
  } else {
    timeinfo = localtime(&time);
  }
  char timecstr[100];
  strftime (timecstr,100,"%Y-%m-%d %H:%M:%S (%Z)",timeinfo);
  //sprintf(timecstr,"%04d-%02d-%02d %02d:%02d:%02d",
  //        timeinfo->tm_year+1900,
  //        timeinfo->tm_mon+1,
  //        timeinfo->tm_mday,
  //        timeinfo->tm_hour,
  //        timeinfo->tm_min,
  //        timeinfo->tm_sec);
  //string timestr = timecstr;
  //if (utc) {
  //  timestr .= " (UTC)";
  //} else {
  //  timestr .= " (local)";
  //}
  return timecstr;
}

/// return a unix-style uname
string CmdLine::unix_uname() const {
  utsname utsbuf;
  int utsret = uname(&utsbuf);
  if (utsret != 0) {return "Error establishing uname";}
  ostringstream uname_result;
  uname_result << utsbuf.sysname << " " 
               << utsbuf.nodename << " "
               << utsbuf.release << " "
               << utsbuf.version << " "
               << utsbuf.machine;
  return uname_result.str();
}

string CmdLine::unix_username() const {
  char * logname;
  logname = getenv("LOGNAME");
  return logname;
}

/// report failure of conversion
void CmdLine::_report_conversion_failure(const string & opt, 
                                         const string & optstring) const {
  ostringstream ostr;
  ostr << "could not convert option ("<<opt<<") value ("
       <<optstring<<") to requested type"<<endl; 
  throw Error(ostr);
}

void CmdLine::assert_all_options_used() const {
  // deal with the help part
  if (__help_enabled) {
    if (present("-h") || present("--help")) {
      print_help();
      exit(0);
    }
  }
  if (! all_options_used()) {
    ostringstream ostr;
    ostr <<"Unrecognised options on the command line" << endl;
    throw Error(ostr);
  }
}


bool CmdLine::Error::_do_printout = true;
CmdLine::Error::Error(const std::ostringstream & ostr) 
  : _message(ostr.str()) {
  if (_do_printout) cerr << "CmdLine Error: " << _message << endl;;
}
CmdLine::Error::Error(const std::string & str) 
  : _message(str) {
  if (_do_printout) cerr << "CmdLine Error: " << _message << endl;;
}

string CmdLine::current_path() const {
  const size_t maxlen = 10000;
  char tmp[maxlen];
  getcwd(tmp,maxlen);
  return string(tmp);
}

string CmdLine::header(const string & prefix) const {
  ostringstream ostr;
  ostr << prefix << "" << command_line() << endl;
  ostr << prefix << "from path: " << current_path() << endl;
  ostr << prefix << "started at: " << time_stamp_at_start() << endl;
  ostr << prefix << "by user: "    << unix_username() << endl;
  ostr << prefix << "running on: " << unix_uname() << endl;
  ostr << prefix << "git state (if any): " << git_info() << endl;
  return ostr.str();
}


string CmdLine::OptionHelp::type_name() const {
  if      (type == typeid(int)   .name()) return "int"   ;
  else if (type == typeid(double).name()) return "double";
  else if (type == typeid(string).name()) return "string";
  else return type;
}

string CmdLine::OptionHelp::summary() const {
  ostringstream ostr;
  if (! required) ostr << "[";
  ostr << option;
  if (takes_value) ostr << " " << argname;
  if (! required) ostr << "]";
  return ostr.str();
}
string CmdLine::OptionHelp::description() const {
  ostringstream ostr;
  ostr << option;
  if (takes_value) {
    ostr << " " << argname << " (" << type_name() << ")";
    if (! required) ostr << "     default: " << default_value;
  }
  ostr << "\n";
  if (help.size() > 0) ostr << "  " << help << endl;
  ostr << endl;
  return ostr.str();
}


void CmdLine::print_help() const {
  // First print a summary
  cout << "\nUsage: \n       " << __arguments[0];
  for (const auto & opt: __options_queried) {
    cout << " " << __options_help[opt].summary();
  }
  cout << endl << endl;

  if (__overall_help_string.size() != 0) {
    cout << __overall_help_string;
    cout << endl << endl;
    cout << "Detailed option help" << endl;
    cout << "--------------------" << endl;
  }
  
  // Then print detailed usage for each option
  for (const auto & opt: __options_queried) {
    const OptionHelp & opthelp = __options_help[opt];
    cout << "  " << opthelp.description();
  }
}

// From https://www.jeremymorgan.com/tutorials/c-programming/how-to-capture-the-output-of-a-linux-command-in-c/
string CmdLine::stdout_from_command(string cmd) const {

  string data;
  FILE * stream;
  const int max_buffer = 1024;
  char buffer[max_buffer];
  cmd.append(" 2>&1");
  
  stream = popen(cmd.c_str(), "r");
  if (stream) {
    while (!feof(stream))
      if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
    pclose(stream);
  }
  return data;
}

//
string CmdLine::git_info() const {
  string log_line = stdout_from_command("git log --pretty='%H %d of %cd' --decorate=short -1");
  for (auto & c : log_line) {if (c == 0x0a || c == 0x0d) c = ';';}
  
  if (log_line.substr(0,6) == "fatal:") {
    log_line = "no git info";
  } else {
    // add info about potentially modified files
    string modifications;
    string status_output = stdout_from_command("git status --porcelain --untracked-files=no");
    for (auto & c : status_output) {if (c == 0x0a || c == 0x0d) c = ',';}
    log_line += "; ";
    log_line += status_output;
  }
  return log_line;
}
