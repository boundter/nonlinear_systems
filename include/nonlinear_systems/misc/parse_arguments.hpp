#ifndef __PARSE_ARGUMENTS__
#define __PARSE_ARGUMENTS__

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <sstream>
#include <fstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;


// https://stackoverflow.com/questions/40754810/default-value-for-vector-valued-boostprogram-options
namespace std
{
  std::ostream& operator<<(std::ostream& os, const std::vector<double>& vec) {
    for (auto item : vec) {
      os << item << " ";
    }
    return os;
  }
  
  
  std::ostream& operator<<(std::ostream& os, const std::vector<unsigned int>& vec) {
    for (auto item : vec) {
      os << item << " ";
    }
    return os;
  }
}


namespace nonlinear_systems {
/*!
 * \brief Base struct for the definition of a vector holding all the arguments.
 */
struct ArgumentBase {
  std::string short_name;
  std::string long_name;
  std::string name;
  std::string description;

  /*!
   *  @param long_name long name of the command line arguments.
   *  @param short_name short name of the command line arguments.
   *  @param description description of the command line argument.
   */
  ArgumentBase(std::string long_name, std::string short_name, std::string description) 
  : short_name(short_name), long_name(long_name), description(description){
    name = long_name + "," + short_name;
  }


  /*!
   *  @param long_name long name of the command line arguments.
   *  @param description description of the command line argument.
   */
  ArgumentBase(std::string long_name, std::string description)
  : long_name(long_name), name(long_name), description(description) {}

  virtual void AddArgument(po::options_description& desc) {}

  virtual void ParseArgument(po::variables_map& vmap) {}

  std::string GetName() {
    return long_name;
  }

  virtual std::string GetValueAsString() {}
};


template<typename T>
/*!
 * \brief Struct to hold the command line arguments.
 */
struct Argument: public ArgumentBase {
  T& value;
  T default_value;

  /*!
   *  @param long_name long name of the command line arguments.
   *  @param description description of the command line argument.
   *  @param value variable holding the argument.
   *  @param default_value default_value of the argument.
   */
  Argument(std::string long_name, std::string description, T& value, T default_value)
  : ArgumentBase(long_name, description), value(value), 
  default_value(default_value) { }


  /*!
   *  @param long_name long name of the command line arguments.
   *  @param short_name short name of the command line arguments.
   *  @param description description of the command line argument.
   *  @param value variable holding the argument.
   *  @param default_value default_value of the argument.
   */
  Argument(std::string long_name, std::string short_name, std::string description, 
      T& value, T default_value)
  : ArgumentBase(long_name, short_name, description), value(value), 
  default_value(default_value) {}


  /*!
   *  Adds the arguments to the description.
   */
  void AddArgument(po::options_description& desc) {
    desc.add_options()
      (name.c_str(), po::value<T>()->default_value(default_value), description.c_str());
  } 


  /*!
   *  Parse the arguments and save the value to the variable. If no value was
   *  passed it will set the default value.
   */
  void ParseArgument(po::variables_map& vmap) {
    if (vmap.count(long_name)) {
      value = vmap[long_name].as<T>();
    }
  }


  std::string GetValueAsString() {
    return std::to_string(value);
  }
};


// specialization for strings
template<> std::string Argument<std::string>
::GetValueAsString() {
  return value;
}


// specialization for vectors
template<> void Argument<std::vector<double>>
::AddArgument(po::options_description& desc) {
  desc.add_options()
    (name.c_str(), po::value<std::vector<double>>()->multitoken()->
     default_value(default_value), description.c_str());
}


template<> void Argument<std::vector<unsigned int>>
::AddArgument(po::options_description& desc) {
  desc.add_options()
    (name.c_str(), po::value<std::vector<unsigned int>>()->multitoken()->
     default_value(default_value), description.c_str());
}


// https://stackoverflow.com/questions/1430757
template<> std::string Argument<std::vector<double>>
::GetValueAsString() {
  std::stringstream ss;
  for (auto it = value.begin(); it != value.end(); ++it) {
    ss << (*it);
    if (it + 1 != value.end()) {
      ss << ",";
    }
  }
  return ss.str();
}


template<> std::string Argument<std::vector<unsigned int>>
::GetValueAsString() {
  std::stringstream ss;
  for (auto it = value.begin(); it != value.end(); ++it) {
    ss << (*it);
    if (it + 1 != value.end()) {
      ss << ",";
    }
  }
  return ss.str();
}


/*!
 *  \brief Parse the command line arguments to variables.
 */
class ParseArguments {
 public:
  std::vector<std::unique_ptr<ArgumentBase>>& arguments;

  /*!
   *  Parse the values from the command line arguments to the variables given as
   *  Arguments. This also build the interface of the help message.
   */
  ParseArguments(int& argc, char* argv[], 
      std::vector<std::unique_ptr<ArgumentBase>>& arguments)
  : arguments(arguments) {
    po::options_description desc("Allowed options");
    AddOptions(desc);
    po::variables_map vmap;
    po::store(po::parse_command_line(argc, argv, desc), vmap);
    po::notify(vmap);
    Parse(desc, vmap);
  }

 protected:
  void AddOptions(po::options_description& desc) {
    desc.add_options()("help,h", "show this help message");
    for (auto it = arguments.begin(); it != arguments.end(); ++it) {
      (*it)->AddArgument(desc);
    }
  }

  
  void Parse(po::options_description& desc, po::variables_map& vmap) {
    if (vmap.count("help")) {
      std::cout << desc;
      exit(0);
    }
    for (auto it = arguments.begin(); it != arguments.end(); ++it) {
      (*it)->ParseArgument(vmap);
    }
  }
};


class WriteArgumentsToFile {
  public:
    std::vector<std::unique_ptr<ArgumentBase>>& arguments;
    FILE* file;

    WriteArgumentsToFile(std::vector<std::unique_ptr<ArgumentBase>>& arguments,
        FILE* outputfile)
      :arguments(arguments), file(outputfile) {
        fprintf(file, "#");
        WriteArguments();
        fprintf(file, "\n");
        fflush(file);
      }

    WriteArgumentsToFile(std::vector<std::unique_ptr<ArgumentBase>>& arguments, 
    std::fstream* file): arguments(arguments) {
      *file << "#";
      WriteArgument(file);
      *file << "\n";
      file->flush();
    }

  protected:

    void WriteArguments() {
      for (auto it = arguments.begin(); it != arguments.end(); ++it) {
        fprintf(file, " %s=%s", (*it)->GetName().c_str(),
            (*it)->GetValueAsString().c_str());
      }
    }

    void WriteArgument(std::fstream* file_stream) {
      for (auto it = arguments.begin(); it != arguments.end(); ++it)
        *file_stream << " " << (*it)->GetName() << "=" << (*it)->GetValueAsString();
    }
};
} // nonlinear_systems
#endif
