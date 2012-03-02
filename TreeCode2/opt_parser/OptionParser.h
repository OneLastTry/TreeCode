/*
 * OptionParser.h
 *
 *  Created on: 2 Mar 2012
 *      Author: stefans
 */

#ifndef OPTIONPARSER_H_
#define OPTIONPARSER_H_

#include <iostream>
#include <boost/lexical_cast.hpp>
#include <string>

class Option{
public:
	Option(std::string longname, std::string shortname, std::string description):
		longname_(longname), shortname_(shortname), description_(description){}

	std::string getLongname() const {return longname_;}
	std::string getShortname() const {return shortname_;}
	std::string getDescription() const {return description_;}
	bool match(char* str){
		return (getLongname().compare(str) == 0 || getShortname().compare(str) == 0);
	}

	void display(std::ostream& out){
		out << longname_ << "," << shortname_ << "\t" << description_ << std::endl;
	}

private:
	std::string longname_, shortname_, description_;
};

class OptionParser{
public:
	OptionParser(int argc, char** argv):argc_(argc), argv_(argv){}

	template <typename T>
	T getOption(Option opt){
		for (int i = 0; i < argc_; i++) {
			if(opt.match(argv_[i])){
				return boost::lexical_cast<T>(argv_[++i]);
			}
		}
		//If no option supplied, print an error and exit
		std::cerr << "Option must be supplied:" << std::endl;
		opt.display(std::cerr);
		exit(1);
	}

	bool optionPresent(Option opt){
		for (int i = 0; i < argc_ ; i++) {
			if(opt.match(argv_[i]))
				return true;
		}
		return false;
	}

private:
	int argc_;
	char** argv_;
};

#endif /* OPTIONPARSER_H_ */
