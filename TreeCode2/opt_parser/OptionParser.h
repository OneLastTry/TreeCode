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
	Option(std::string longname, std::string shortname, std::string description, bool compulsory):
		longname_(longname), shortname_(shortname), description_(description), compulsory_(compulsory){}

	std::string getLongname() const {return longname_;}
	std::string getShortname() const {return shortname_;}
	std::string getDescription() const {return description_;}
	bool match(char* str) const{
		return (getLongname().compare(str) == 0 || getShortname().compare(str) == 0);
	}

	void display(std::ostream& out) const{
		out << longname_ << "," << shortname_ << "\t" << description_ << std::endl;
	}

	bool isCompulsory() const{
		return compulsory_;
	}

	virtual void set(char** argv, int& i) = 0;

private:
	std::string longname_, shortname_, description_;
	bool compulsory_;
};

template <typename T>
class ArgOption : public Option {
public:
	ArgOption(std::string longname, std::string shortname, std::string description, T& reference, bool compulsory = false):
		Option(longname, shortname, description, compulsory), reference_(reference){}

	void set(char** argv, int& i){
		reference_ = boost::lexical_cast<T>(argv[++i]);
	}
private:
	T& reference_;
};

class BoolOption : public Option {
public:
	BoolOption(std::string longname, std::string shortname, std::string description, bool& reference, bool compulsory = false):
		Option(longname, shortname, description, compulsory), reference_(reference){
		reference_ = false;
	}

	void set(char** argv, int& i){
		reference_ = true;
	}
private:
	bool& reference_;
};

class OptionParser{
public:
	OptionParser(int argc, char** argv):argc_(argc), argv_(argv){}

	~OptionParser(){
		for(std::vector<Option*>::iterator it = options.begin(); it < options.end(); it++)
			delete (*it);
	}

	OptionParser& operator<<(Option* o){
		options.push_back(o);
		return *this;
	}

	void display(std::ostream& out = std::cout){
		for(std::vector<Option*>::iterator it = options.begin(); it < options.end(); it++)
			(*it)->display(out);
	}

	unsigned int parse(){
		unsigned int num_ops = 0;

		for(std::vector<Option*>::iterator it = options.begin(); it < options.end(); it++){
			bool found = false;
			for(int i=0;i<argc_;i++){
				if((*it)->match(argv_[i])){
					(*it)->set(argv_, i);
					num_ops++;
					found = true;
					break;
				}
			}
			//If we reached here without finding it, then check if it /must/ exist
			if((*it)->isCompulsory() && !found){
				std::cerr << "You must set this option:" << std::endl;
				(*it)->display(std::cerr);
				std::cerr << "Usage:" << std::endl;
				display(std::cerr);
				exit(1);
			}
		}
		return num_ops;
	}

	unsigned int size() const {return options.size();}

private:
	int argc_;
	char** argv_;
	std::vector<Option*> options;
};

#endif /* OPTIONPARSER_H_ */
