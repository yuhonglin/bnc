/*
  An adapter to extend std::string
 */

#ifndef STRING_H
#define STRING_H

#include <string>

#include <boost/algorithm/string.hpp>

namespace bnc {

    class String
    {
	string s_;
    public:
	String(const std::string& s) s_(s) {};

	String upper() const {
	    return boost::to_upper_copy(s_);
	}

	String lower() const {
	    return boost::to_lower_copy(s_);
	}
	
	vector<String> split(const string& delimiter) const {
	    vector<string> parts;
	    boost::split(parts, s_, boost::is_any_of(delimiter));
	    return vector<String>(parts.begin(), parts.end());
	}

	String lstrip(const std::string &s) const {
	    boost::trim_left_copy(s_, boost::is_any_of(s));
	}

	String rstrip(const std::string &s) const {
	    boost::trim_right_copy(s_, boost::is_any_of(s));
	}

	String strip(const std::string &s) const {
	    boost::trim_copy(s_, boost::is_any_of(s));
	}

	size_t len() const { return s_.length(); }
	
	~String() {};
    };

    
}  // namespace bnc

#endif /* STRING_H */
