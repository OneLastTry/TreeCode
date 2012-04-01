/*
 * OptionParser.h
 *
 *  Created on: 31 Mar 2012
 *      Author: stefans
 */

#ifndef OPTIONPARSER_H_
#define OPTIONPARSER_H_

#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

namespace po = boost::program_options;

/**
 * @brief Called to read into an Eigen::VectorXd from the format '(x,y,z)'
 * @param in		Input stream.
 * @param[out] vec 	Output vector
 * @return @p in
 */
std::istream& operator>>(std::istream& in, Eigen::VectorXd& vec){
	po::validation_error err = po::validation_error(po::validation_error::invalid_option_value);

	//Set size to zero initially
    vec.resize(0);
    double d;
    int index = 0;

    //Make sure it starts with '('
    int c = in.get();
    if(c != '(')
    	throw err;
    while(true){
    	//Read double, add to vector
        in >> d;
        vec.conservativeResize(index + 1);
        vec[index++] = d;

        //Require ')' or ',' -- ')' ends
        //Would read whitespace, but boost doesn't allow it anyway.
        c = in.get();
        if(c == ')')break;
        else if(c != ',')
        	throw err;
    }
    return in;
}

/**
 * @brief Option parser, wrapping boost's program_options.
 *
 * This class extends boost::program_options::options_description,
 * so the usual methods apply. Also, it implements methods for
 * reading the values back, and a custom parser for Eigen::VectorXd
 * objects.
 *
 * If setUsesConfigFile() is called, then a "--config-file" option
 * is added to the command line arguments in a "config" section.
 * Then, when parse() is called, OptionParser first reads the name
 * of the config file.
 *
 * parse() then reads the command line options, followed by the config
 * file options, so that the command line options get priority.
 */
class OptionParser : public po::options_description{
public:
	/**
	 * By default, add a help option.
	 * @param caption
	 */
	OptionParser(std::string caption, bool addhelp = true) :
		po::options_description(caption), use_config_(false){
		if(addhelp)
			this->add_options()("help,h", "Display this help text");
	}

	~OptionParser(){}

	void parse(int argc, char **argv){
		try{
			//Read the command line options
			po::store(po::parse_command_line(argc, argv, *this), vm_);
			if(vm_.count("help")){
				std::cout << *this << std::endl;
				exit(0);
			}

			//If we are using a config file as well, read from that
			if(use_config_){
				//Get config file name
				if(vm_.count("config-file")){
					std::string config_file_name = vm_["config-file"].as<std::string>();

					//Read from config file
					std::ifstream config_istream(config_file_name.c_str());
					po::store(po::parse_config_file(config_istream, *this), vm_);
				}
			}
			//Call notify -- this will cause exceptions to be thrown, etc
			po::notify(vm_);
		}catch(po::required_option& re){
			std::cerr << "Required option '" << re.get_option_name() << "' not supplied." << std::endl;
			std::cerr << "Use --help for a list of all options." << std::endl;
			exit(1);
		}catch(po::invalid_option_value& ie){
			std::cerr << ie.what() << std::endl;
			exit(1);
		}
		catch(po::error& e){
			std::cerr << e.what() << std::endl;
			exit(1);
		}
	}

	template <typename T>
	const T& get(const std::string& name) const{
		return vm_[name].as<T>();
	}

	void setUsesConfigFile(bool use_config){
		use_config_ = use_config;
		if(use_config){
			//Add an option, under a "Configuration" section,
			//that allows the filename to be specified.
			po::options_description config_opts("Configuration Options");
			config_opts.add_options()
					("config-file", po::value<std::string>(), "Configuration file to read options from");
			this->add(config_opts);
		}
	}

	int count(std::string s) const {
		return vm_.count(s);
	}
private:
	bool use_config_;
	po::variables_map vm_;
};

#endif /* OPTIONPARSER_H_ */
