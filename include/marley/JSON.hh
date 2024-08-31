/// @file
/// @copyright Copyright (C) 2016-2024 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

// Based on https://github.com/nbsdx/SimpleJSON

#pragma once

// standard library includes
#include <cctype>
#include <cmath>
#include <cstdint>
#include <deque>
#include <initializer_list>
#include <iostream>
#include <istream>
#include <ostream>
#include <limits>
#include <sstream>
#include <map>
#include <string>
#include <type_traits>

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/FileManager.hh"
#include "marley/Logger.hh"

namespace marley {

  namespace {
    std::string json_escape( const std::string& str ) {
      std::string output;
      for( unsigned i = 0; i < str.length(); ++i )
      switch( str[i] ) {
        case '\"': output += "\\\""; break;
        case '\\': output += "\\\\"; break;
        case '\b': output += "\\b";  break;
        case '\f': output += "\\f";  break;
        case '\n': output += "\\n";  break;
        case '\r': output += "\\r";  break;
        case '\t': output += "\\t";  break;
        default  : output += str[i]; break;
      }
      return output;
    }
  }

  class JSON
  {
    union Data {
      Data(double d) : float_(d) {}
      Data(long l) : integer_(l) {}
      Data(bool b) : boolean_(b) {}
      Data(const std::string& s) : string_(new std::string(s)) {}
      Data() : integer_(0) {}

      std::deque<JSON>* list_;
      std::map<std::string, JSON>* map_;
      std::string* string_;
      double float_;
      long integer_;
      bool boolean_;
    } data_;

    public:
      enum class DataType {
        Null,
        Object,
        Array,
        String,
        Floating,
        Integral,
        Boolean
      };

      template <typename Container> class JSONWrapper {

        private:
          Container* object;

        public:
          JSONWrapper(Container* val) : object(val) {}
          JSONWrapper(std::nullptr_t) : object(nullptr) {}

          typename Container::iterator begin() {
            return object ? object->begin() : typename Container::iterator();
          }
          typename Container::iterator end() {
            return object ? object->end() : typename Container::iterator();
          }
          typename Container::const_iterator begin() const {
            return object ? object->begin() : typename Container::iterator();
          }
          typename Container::const_iterator end() const {
            return object ? object->end() : typename Container::iterator();
          }
      };

      template <typename Container>
      class JSONConstWrapper {

        private:
          const Container* object;

        public:
          JSONConstWrapper(const Container* val) : object(val) {}
          JSONConstWrapper(std::nullptr_t) : object(nullptr) {}

          typename Container::const_iterator begin() const
            { return object ? object->begin()
              : typename Container::const_iterator(); }
          typename Container::const_iterator end() const
            { return object ? object->end()
              : typename Container::const_iterator(); }
      };

      JSON() : data_(), type_(DataType::Null) {}

      JSON(std::initializer_list<JSON> list) : JSON()
      {
        set_type(DataType::Object);
        for(auto i = list.begin(), e = list.end(); i != e; ++i, ++i)
          operator[](i->to_string()) = *std::next(i);
      }

      JSON(JSON&& other) : data_(other.data_), type_(other.type_)
      { other.type_ = DataType::Null; other.data_.map_ = nullptr; }

      JSON& operator=(JSON&& other) {
        data_ = other.data_;
        type_ = other.type_;
        other.data_.map_ = nullptr;
        other.type_ = DataType::Null;
        return *this;
      }

      JSON(const JSON& other) {
        switch(other.type_) {
        case DataType::Object:
          data_.map_ = new std::map<std::string,JSON>(
            other.data_.map_->begin(), other.data_.map_->end());
          break;
        case DataType::Array:
          data_.list_ = new std::deque<JSON>(other.data_.list_->begin(),
            other.data_.list_->end());
          break;
        case DataType::String:
          data_.string_ = new std::string(*other.data_.string_);
          break;
        default:
          data_ = other.data_;
        }
        type_ = other.type_;
      }

      JSON& operator=(const JSON& other) {
        switch(other.type_) {
          case DataType::Object:
            data_.map_ = new std::map<std::string,JSON>(
              other.data_.map_->begin(), other.data_.map_->end());
            break;
          case DataType::Array:
            data_.list_ = new std::deque<JSON>( other.data_.list_->begin(),
              other.data_.list_->end());
            break;
          case DataType::String:
            data_.string_ = new std::string(*other.data_.string_);
            break;
          default:
            data_ = other.data_;
        }
        type_ = other.type_;
        return *this;
      }

      ~JSON() {
        switch(type_) {
          case DataType::Array:
            delete data_.list_;
            break;
          case DataType::Object:
            delete data_.map_;
            break;
          case DataType::String:
            delete data_.string_;
            break;
          default:;
        }
      }

      template <typename T> JSON(T b,
        typename std::enable_if<std::is_same<T,bool>::value>::type* = 0)
        : data_(b), type_(DataType::Boolean) {}

      template <typename T> JSON(T i,
        typename std::enable_if<std::is_integral<T>::value
        && !std::is_same<T,bool>::value>::type* = 0)
        : data_(static_cast<long>(i)), type_(DataType::Integral) {}

      template <typename T> JSON(T f,
        typename std::enable_if<std::is_floating_point<T>::value>::type* = 0)
        : data_(static_cast<double>(f)), type_(DataType::Floating) {}

      explicit JSON(const std::string& s)
        : data_(s), type_(DataType::String) {}

      //template <typename T> JSON(T s,
      //  typename std::enable_if<std::is_convertible<T,
      //  std::string>::value>::type* = 0) : data_(std::string(s)),
      //  type_(DataType::String) {}

      JSON(std::nullptr_t) : data_(), type_(DataType::Null) {}

      static inline JSON make(DataType type) {
        JSON ret;
        ret.set_type(type);
        return ret;
      }

      static inline JSON array() {
        return JSON::make(JSON::DataType::Array);
      }

      template <typename... T>
        static JSON array( T... args )
      {
        JSON arr = JSON::make(JSON::DataType::Array);
        arr.append(args...);
        return arr;
      }

      static inline JSON object() {
        return JSON::make(JSON::DataType::Object);
      }

      static inline JSON load(const std::string& s);
      static inline JSON load(std::istream& is);
      static inline JSON load_file(const std::string& s);

      template <typename T> void append(T arg) {
        set_type(DataType::Array);
        data_.list_->emplace_back(arg);
      }

      template <typename T, typename... U> void append(T arg, U... args) {
        append(arg); append(args...);
      }

      template <typename T>
        typename std::enable_if<std::is_same<T,bool>::value, JSON&>::type
        operator=(T b)
      {
          set_type(DataType::Boolean);
          data_.boolean_ = b;
          return *this;
      }

      template <typename T> typename std::enable_if<std::is_integral<T>::value
        && !std::is_same<T,bool>::value, JSON&>::type operator=(T i)
      {
        set_type( DataType::Integral );
        data_.integer_ = i;
        return *this;
      }

      template <typename T>
        typename std::enable_if<std::is_floating_point<T>::value, JSON&>::type
        operator=(T f)
      {
        set_type(DataType::Floating);
        data_.float_ = f;
        return *this;
      }

      template <typename T> typename std::enable_if<std::is_convertible<T,
        std::string>::value, JSON&>::type operator=(T s)
      {
        set_type(DataType::String);
        *data_.string_ = std::string(s);
        return *this;
      }

      JSON& operator[](const std::string& key) {
        set_type(DataType::Object);
        return data_.map_->operator[](key);
      }

      JSON& operator[](unsigned index) {
        set_type(DataType::Array);
        if (index >= data_.list_->size()) data_.list_->resize(index + 1);
        return data_.list_->operator[](index);
      }

      JSON& at(const std::string& key) {
        return operator[](key);
      }

      const JSON& at(const std::string &key) const {
        return data_.map_->at(key);
      }

      JSON& at(unsigned index) {
        return operator[](index);
      }

      const JSON& at(unsigned index) const {
        return data_.list_->at(index);
      }

      int length() const {
        if (type_ == DataType::Array) return data_.list_->size();
        else return -1;
      }

      bool has_key(const std::string& key) const {
        if (type_ == DataType::Object)
          return data_.map_->find( key ) != data_.map_->end();
        else return false;
      }

      int size() const {
        if (type_ == DataType::Object)
          return data_.map_->size();
        else if (type_ == DataType::Array)
          return data_.list_->size();
        else
          return -1;
      }

      inline DataType type() const { return type_; }

      /// Functions for getting primitives from the JSON object.
      inline bool is_null() const { return type_ == DataType::Null; }
      inline bool is_object() const { return type_ == DataType::Object; }
      inline bool is_array() const { return type_ == DataType::Array; }
      inline bool is_string() const { return type_ == DataType::String; }
      inline bool is_float() const { return type_ == DataType::Floating; }
      inline bool is_integer() const { return type_ == DataType::Integral; }
      inline bool is_bool() const { return type_ == DataType::Boolean; }

      std::string to_string() const {
        bool b;
        return to_string(b);
      }

      std::string to_string(bool& ok) const {
        ok = (type_ == DataType::String);
        return ok ? json_escape(*data_.string_) : std::string("");
      }

      std::string to_string_or_throw() const {
        bool ok;
        std::string result = to_string(ok);
        if (!ok) throw marley::Error("Failed to convert JSON value to string");
        return result;
      }

      double to_double() const {
        bool b;
        return to_double(b);
      }

      double to_double(bool& ok) const {
        ok = (type_ == DataType::Floating);
        if (ok) return data_.float_;
        ok = (type_ == DataType::Integral);
        if (ok) return data_.integer_;
        return 0.;
      }

      double to_double_or_throw() const {
        bool ok;
        double result = to_double(ok);
        if (!ok) throw marley::Error("Failed to convert JSON value '"
          + to_string() + "' to double");
        return result;
      }

      long to_long() const {
        bool b;
        return to_long( b );
      }

      long to_long(bool& ok) const {
        ok = (type_ == DataType::Integral);
        return ok ? data_.integer_ : 0;
      }

      long to_long_or_throw() const {
        bool ok;
        double result = to_long(ok);
        if (!ok) throw marley::Error("Failed to convert JSON value '"
          + to_string() + "' to long");
        return result;
      }

      bool to_bool() const {
        bool b;
        return to_bool( b );
      }

      bool to_bool(bool& ok) const {
        ok = (type_ == DataType::Boolean);
        return ok ? data_.boolean_ : false;
      }

      bool to_bool_or_throw() const {
        bool ok;
        double result = to_bool(ok);
        if (!ok) throw marley::Error("Failed to convert JSON value '"
          + to_string() + "' to bool");
        return result;
      }

      JSONWrapper<std::map<std::string,JSON> > object_range() {
        if (type_ == DataType::Object)
          return JSONWrapper<std::map<std::string,JSON>>(data_.map_);
        else return JSONWrapper<std::map<std::string,JSON>>(nullptr);
      }

      JSONWrapper<std::deque<JSON> > array_range() {
        if (type_ == DataType::Array)
          return JSONWrapper<std::deque<JSON>>(data_.list_);
        else return JSONWrapper<std::deque<JSON>>(nullptr);
      }

      JSONConstWrapper<std::map<std::string,JSON> > object_range() const {
        if (type_ == DataType::Object)
          return JSONConstWrapper<std::map<std::string,JSON>>(data_.map_);
        else return JSONConstWrapper<std::map<std::string,JSON>>(nullptr);
      }


      JSONConstWrapper<std::deque<JSON>> array_range() const {
        if ( type_ == DataType::Array )
          return JSONConstWrapper<std::deque<JSON>>(data_.list_);
        else return JSONConstWrapper<std::deque<JSON>>(nullptr);
      }

      // Portions of the serialization functions (dump_string, print)
      // are based on techniques used in the JSON for Modern C++
      // library by Niels Lohmann (https://github.com/nlohmann/json).
      std::string dump_string(const int indent_step = -1) const {
        std::stringstream out;
        // Enable pretty-printing if the user specified a nonnegative
        // indent_step value
        if (indent_step >= 0)
          print(out, static_cast<unsigned int>(indent_step), true);
        // Otherwise, print the JSON object in the most compact form possible
        else print(out, 0, false);

        // Return the completed JSON string
        return out.str();
      }

      // Implementation of serialization to text. Used by the public
      // dump_string() method.
      void print(std::ostream& out, const unsigned int indent_step,
        bool pretty, const unsigned int current_indent = 0) const
      {
        // Use max_digits10 for outputting double-precision floating-point
        // numbers. This ensures that repeated input/output via JSON will
        // not result in any loss of precision. For more information, please
        // see http://tinyurl.com/p8wyhnn
        static std::ostringstream out_float;
        static bool set_precision = false;
        if (!set_precision) {
          out_float.precision(std::numeric_limits<double>::max_digits10);
          set_precision = true;
        }

        unsigned int indent = current_indent;

        switch( type_ ) {
          case DataType::Null:
            out << "null";
            return;
          case DataType::Object: {
            out << '{';
            if (pretty) {
              indent += indent_step;
              out << '\n';
            }
            bool skip = true;
            for( auto &p : *data_.map_ ) {
              if ( !skip ) {
                out << ',';
                if (pretty) out << '\n';
              }

              out << std::string(indent, ' ') << '\"'
                << json_escape( p.first ) << '\"';

              if (pretty) out << " : ";
              else out << ':';

              p.second.print( out, indent_step, pretty, indent );
              skip = false;
            }
            if (pretty) {
              indent -= indent_step;
              out << '\n';
            }
            out << std::string(indent, ' ') + '}';
            return;
          }
          case DataType::Array: {
            out << '[';
            if (pretty) {
              indent += indent_step;
              out << '\n';
            }
            bool skip = true;
            for( auto &p : *data_.list_ ) {
              if ( !skip ) {
                out << ',';
                if (pretty) out << '\n';
              }
              out << std::string(indent, ' ');
              p.print( out, indent_step, pretty, indent );
              skip = false;
            }
            if (pretty) {
              indent -= indent_step;
              out << '\n';
            }
            out << std::string(indent, ' ') << ']';
            return;
          }
          case DataType::String:
            out << '\"' + json_escape( *data_.string_ ) + '\"';
            return;
          case DataType::Floating:
            // Clear any previous contents of the stringstream
            out_float.str("");
            out_float.clear();
            // Fill it with the new floating-point number
            out_float << data_.float_;
            // Output the resulting string to the stream
            out << out_float.str();
            return;
          case DataType::Integral:
            out << data_.integer_;
            return;
          case DataType::Boolean:
            out << (data_.boolean_ ? "true" : "false");
            return;
          default:
            break;
        }

        return;
      }

      void check_if_object( const std::string& key ) const {
        if ( type_ != DataType::Object ) throw marley::Error( "Attempted"
          " to retrieve a value for the key '" + key + "' from a JSON"
          " primitive that is not an object" );
      }

    private:

      void set_type( DataType type ) {
        if ( type == type_ ) return;

        switch( type_ ) {
          case DataType::Object:
            delete data_.map_;
            break;
          case DataType::Array:
            delete data_.list_;
            break;
          case DataType::String:
            delete data_.string_;
            break;
          default:;
        }

        switch( type ) {
          case DataType::Null:
            data_.map_ = nullptr;
            break;
          case DataType::Object:
            data_.map_ = new std::map< std::string, JSON >();
            break;
          case DataType::Array:
            data_.list_ = new std::deque< JSON >();
            break;
          case DataType::String:
            data_.string_ = new std::string();
            break;
          case DataType::Floating:
            data_.float_ = 0.;
            break;
          case DataType::Integral:
            data_.integer_ = 0;
            break;
          case DataType::Boolean:
            data_.boolean_ = false;
            break;
        }

        type_ = type;
      }

    public:

      // Attempts to get a floating point number from a JSON object with
      // a given key. If the attempt fails, throw a marley::Error.
      double get_double( const std::string& key ) const {
        check_if_object( key );
        if ( has_key(key) ) return this->at( key ).to_double_or_throw();
        else throw marley::Error( "Missing JSON key '" + key + '\'' );
        return 0.;
      }

      // Attempts to get a floating point number from a JSON object with
      // a given key. If the key doesn't exist, use a default value. If a
      // conversion attempt fails, throw a marley::Error.
      double get_double( const std::string& key, double default_value ) const {
        check_if_object( key );
        if ( !has_key(key) ) return default_value;
        else return this->at( key ).to_double_or_throw();
      }

      // Attempts to get an integer from a JSON object with
      // a given key. If the attempt fails, throw a marley::Error.
      long get_long( const std::string& key ) const {
        check_if_object( key );
        if ( has_key(key) ) return this->at( key ).to_long_or_throw();
        else throw marley::Error( "Missing JSON key '" + key + '\'' );
        return 0.;
      }

      // Attempts to get an integer from a JSON object with
      // a given key. If the key doesn't exist, use a default value. If a
      // conversion attempt fails, throw a marley::Error.
      long get_long( const std::string& key, long default_value ) const {
        check_if_object( key );
        if ( !has_key(key) ) return default_value;
        else return this->at( key ).to_long_or_throw();
      }

      // Attempts to get a bool from a JSON object with
      // a given key. If the attempt fails, throw a marley::Error.
      bool get_bool( const std::string& key ) const {
        check_if_object( key );
        if ( has_key(key) ) return this->at( key ).to_bool_or_throw();
        else throw marley::Error( "Missing JSON key '" + key + '\'' );
        return 0.;
      }

      // Attempts to get a bool from a JSON object with
      // a given key. If the key doesn't exist, use a default value. If a
      // conversion attempt fails, throw a marley::Error.
      bool get_bool( const std::string& key, bool default_value ) const {
        check_if_object( key );
        if ( !has_key(key) ) return default_value;
        else return this->at( key ).to_bool_or_throw();
      }

      // Attempts to get a string from a JSON object with
      // a given key. If the attempt fails, throw a marley::Error.
      std::string get_string( const std::string& key ) const {
        check_if_object( key );
        if ( has_key(key) ) return this->at( key ).to_string_or_throw();
        else throw marley::Error( "Missing JSON key '" + key + '\'' );
        return std::string( "" );
      }

      // Attempts to get a string from a JSON object with
      // a given key. If the key doesn't exist, use a default value. If a
      // conversion attempt fails, throw a marley::Error.
      std::string get_string( const std::string& key,
        const std::string& default_value ) const
      {
        check_if_object( key );
        if ( !has_key(key) ) return default_value;
        else return this->at( key ).to_string_or_throw();
      }

      // Copies a subobject from a JSON object with a given key. If the attempt
      // fails, throw a marley::Error, unless the user asks us not to do so.
      marley::JSON get_object( const std::string& key,
        bool throw_error = true ) const
      {
        check_if_object( key );
        if ( has_key(key) ) return this->at( key );
        else if ( throw_error ) marley::Error(
          "Missing JSON key '" + key + '\'' );
        return JSON::make( JSON::DataType::Object );
      }

    private:

      DataType type_ = DataType::Null;
  };

  namespace {

    JSON parse_next( std::istream& );

    void issue_parse_error( char found_char, const std::string& message )
    {
      std::string msg( message );
      if ( found_char == std::ifstream::traits_type::eof() )
        msg += "end-of-file";
      else msg += std::string( "\'" ) + found_char + '\'';
      //MARLEY_LOG_WARNING() << msg;
      throw marley::Error( msg );
    }

    void issue_parse_error( const std::string& found_str,
      const std::string& message, const std::istream& is )
    {
      std::string msg( message );
      if ( !is ) msg += "end-of-file";
      else msg += '\'' + found_str + '\'';
      //MARLEY_LOG_WARNING() << msg;
      throw marley::Error( msg );
    }

    // Skips single-line comments // and multi-line comments /* */
    // These are technically not valid in JSON (the standard doesn't allow
    // comments), but they are valid in Javascript object literals.
    void skip_comment( std::istream& in, bool is_multiline = false ) {
      if ( is_multiline ) {
        char c;
        while ( in.get(c) ) {
          if ( c == '*' && in.peek() == '/' ) {
            in.ignore();
            break;
          }
        }
      }
      // Ignore all further characters until either a newline or end-of-file
      else in.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    }

    // Skips whitespace and comments, saving the last character read to
    // read_char.
    void skip_ws( std::istream& in, char& read_char ) {
      while ( read_char = in.get(), std::isspace(read_char) ) continue;
      if ( read_char == '/' ) {
        char c = in.peek();
        if ( c == '/' || c == '*' ) {
          read_char = in.get();
          skip_comment( in, c == '*' );
          return skip_ws( in, read_char );
        }
      }
    }

    // Removes whitespace and comments from the input stream, putting back
    // the first non-whitespace and non-comment character it finds.
    void consume_ws( std::istream& in ) {
      static char next;
      skip_ws( in, next );
      in.putback( next );
    }

    // Removes whitespace and comments from the input stream, returning the
    // first non-whitespace and non-comment character it finds.
    char get_next_char( std::istream& in )
    {
      static char next;
      skip_ws( in, next );
      return next;
    }

    JSON parse_object( std::istream& in ) {

      JSON object = JSON::make( JSON::DataType::Object );

      for ( ;; ) {

        consume_ws( in );
        JSON key;

        if ( in.peek() == '}' ) {
          in.ignore();
          return object;
        }
        else if ( in.peek() == '\"' ) {
          key = parse_next(in);
        }
        // The key isn't quoted, so assume it's a single word followed
        // by a colon. Note that vanilla JSON requires all keys to be quoted,
        // but Javascript object literals allow unquoted keys.
        else {
          std::string key_str;
          char c;
          while ( in.get(c) ) {
            if ( c == ':' || std::isspace(c) ) {
              in.putback( c );
              break;
            }
            key_str += c;
          }
          key = key_str;
        }

        char next = get_next_char( in );
        if ( next != ':' ) {
          issue_parse_error( next, "JSON object: Expected colon, found " );
          break;
        }

        consume_ws( in );
        JSON value = parse_next( in );
        object[ key.to_string() ] = value;

        next = get_next_char( in );
        if ( next == ',' ) continue;
        else if ( next == '}' ) break;
        else {
          issue_parse_error( next, "JSON object: Expected comma, found " );
          break;
        }
      }

      return object;
    }

    JSON parse_array(std::istream& in) {
      JSON array = JSON::make(JSON::DataType::Array);
      unsigned index = 0;

      for (;;) {

        consume_ws(in);
        if (in.peek() == ']') {
          in.ignore();
          return array;
        }

        array[index++] = parse_next(in);
        consume_ws(in);

        char next = in.get();
        if (next == ',') continue;
        else if (next == ']') break;
        else {
          issue_parse_error(next, "JSON array: Expected ',' or ']'"
            ", found ");
          return JSON::make(JSON::DataType::Array);
        }
      }

      return array;
    }

    JSON parse_string(std::istream& in) {
      JSON str;
      std::string val;
      for(char c = in.get(); c != '\"' && in; c = in.get()) {
        if (c == '\\') {
          switch( in.get() ) {
            case '\"': val += '\"'; break;
            case '\\': val += '\\'; break;
            case '/' : val += '/' ; break;
            case 'b' : val += '\b'; break;
            case 'f' : val += '\f'; break;
            case 'n' : val += '\n'; break;
            case 'r' : val += '\r'; break;
            case 't' : val += '\t'; break;
            case 'u' : {
              val += "\\u" ;
              for(unsigned i = 1; i <= 4; ++i) {
                c = in.get();
                if ((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')
                  || (c >= 'A' && c <= 'F')) val += c;
                else {
                  issue_parse_error(c, "JSON string: Expected hex character"
                    " in unicode escape, found ");
                  return JSON::make(JSON::DataType::String);
                }
              }
              break;
            }
            default: val += '\\'; break;
          }
        }
        else val += c;
      }
      str = val;
      return str;
    }

    //FIXME
    JSON parse_number(std::istream& in, char old) {
      JSON Number;
       std::string val, exp_str;
      char c = old;
      bool isDouble = false;
      long exp = 0;
      for (;;) {
        if ( (c == '-') || (c >= '0' && c <= '9') )
          val += c;
        else if ( c == '.' ) {
          val += c;
          isDouble = true;
        }
        else
          break;
        c = in.get();
      }
      if ( c == 'E' || c == 'e' ) {
        if ( in.peek() == '-' ) { in.ignore(); exp_str += '-'; }
        else if ( in.peek() == '+' ) { in.ignore(); }
        for (;;) {
          c = in.get();
          if ( c >= '0' && c <= '9' )
            exp_str += c;
          else if ( !std::isspace( c ) && c != ',' && c != ']' && c != '}' ) {
            issue_parse_error(c, "JSON number: Expected a number for"
              " exponent, found ");
            return JSON::make(JSON::DataType::Null);
          }
          else
            break;
        }
        exp = std::stol( exp_str );
      }
      else if ( !std::isspace( c ) && c != ',' && c != ']' && c != '}' ) {
        issue_parse_error(c, "JSON number: unexpected character ");
        return JSON::make(JSON::DataType::Null);
      }
      in.putback(c);

      if ( isDouble )
        Number = std::stod( val ) * std::pow( 10, exp );
      else {
        if ( !exp_str.empty() )
          Number = std::stol( val ) * std::pow( 10, exp );
        else
          Number = std::stol( val );
      }
      return  Number ;
    }

    JSON parse_bool(std::istream& in, char old) {
      JSON b;
      std::string s(1, old);
      if (old == 't') {
        for (size_t i = 0; i < 3; ++i) s += in.get();
        if (s == "true") b = true;
      }
      else if (old == 'f') {
        for (size_t i = 0; i < 4; ++i) s += in.get();
        if (s == "false") b = false;
      }
      if (b.type() == JSON::DataType::Null) {
        // Get the entire string if the user supplied an invalid value
        while (in.good() && !std::isspace(in.peek())) s += in.get();
        marley_utils::trim_inplace(s);

        issue_parse_error(s, "JSON bool: Expected 'true' or 'false', found ",
          in);
        return JSON::make(JSON::DataType::Null);
      }
      return b;
    }

    JSON parse_null(std::istream& in) {
      JSON null;
      std::string s(1, 'n');
      for (size_t i = 0; i < 3; ++i) s += in.get();
      if ( s != "null") {
        issue_parse_error("JSON null: Expected 'null', found ", s, in);
        return JSON::make(JSON::DataType::Null);
      }
      return null;
    }

    JSON parse_include( std::istream& in ) {
      std::string s( 1, '#' );
      for (size_t i = 0; i < 9; ++i) s += in.get();
      if ( s != "#include:\"") {
        throw marley::Error( "JSON include: Expected 'include:\"', found '"
          + s + '\'' );
        return JSON::make( JSON::DataType::Null );
      }

      // Parse the included file name into a temporary JSON object, then find
      // the full path to the file
      JSON file_name_json = parse_string( in );
      std::string file_name = file_name_json.to_string();

      const auto& fm = marley::FileManager::Instance();
      std::string full_file_name = fm.find_file( file_name );

      if ( full_file_name.empty() ) {
        throw marley::Error( "Could not locate the included JSON file \""
          + file_name + "\". Please check that the file name is spelled"
          " correctly and that the file is in a folder on the MARLEY"
          " search path." );
      }

      // Open the file for reading and check that it is ready to use
      std::ifstream included_file_stream( full_file_name );
      if ( !included_file_stream.good() ) {
        throw marley::Error( "Could not read from the included JSON file \""
          + full_file_name + '\"' );
      }

      // Use a recursive call to parse_next() to interpret the JSON in the
      // file, allowing for the possibility of nested #include commands
      return parse_next( included_file_stream );
    }

    JSON parse_next( std::istream& in ) {
      char value = get_next_char( in );
      switch(value) {
        case '[' : return parse_array(in);
        case '{' : return parse_object(in);
        case '\"': return parse_string(in);
        case 't' :
        case 'f' : return parse_bool(in, value);
        case 'n' : return parse_null(in);
        case '#' : return parse_include(in);
        default  :
          if ((value <= '9' && value >= '0') || value == '-')
            return parse_number(in, value);
      }
      // Complain and throw an error if there was a problem
      if (!in) throw marley::Error("Unexpected end of JSON configuration"
        " file found\n");
      else throw marley::Error(std::string("JSON parse:")
        + " Unknown starting character '" + value + "'\n");
      return JSON();
    }
  }

  inline JSON JSON::load_file(const std::string& filename) {
    std::ifstream in(filename);
    if (in.good()) return load(in);
    else {
      throw marley::Error("Could not open the file \"" + filename + "\"");
      return JSON::make(JSON::DataType::Null);
    }
  }

  inline JSON JSON::load(std::istream& in) {
    char first = get_next_char( in );
    if (first != '{') {
      throw marley::Error("Missing '{' at beginning of JSON object");
      in.putback(first);
      return parse_object(in);
    }
    else {
      in.putback(first);
      return parse_next(in);
    }
  }

  inline JSON JSON::load(const std::string& str) {
    std::stringstream iss(str);
    return load(iss);
  }

}

// Stream operators for JSON input and output using C++ streams
inline std::ostream& operator<<(std::ostream &os, const marley::JSON& json) {
  os << json.dump_string();
  return os;
}

inline std::istream& operator>>(std::istream &is, marley::JSON& json) {
  json = marley::JSON::load(is);
  return is;
}

/////////////////////////////////////////////////////////
// Utility function templates for retrieval of JSON data
/////////////////////////////////////////////////////////

// Detects whether type T has a member function 'push_back' that is callable
// with an argument of T::value_type. Examples include std::vector,
// std::deque, etc.
template < typename T, typename = void > struct HasPushBack
  : std::false_type {};

template < typename T > struct HasPushBack< T, std::void_t< decltype(
  std::declval< T >().push_back( std::declval< typename T::value_type >() )
  ) > > : std::true_type {};

// Detects whether type T has a member function 'clear' that is callable with
// an argument of T::value_type. Examples include std::vector, std::deque,
// etc.
template < typename T, typename = void > struct HasClear
  : std::false_type {};

template < typename T > struct HasClear< T, std::void_t< decltype(
  std::declval< T >().clear() ) > > : std::true_type {};

// Detects whether type T has member types 'key_type' and 'mapped_type' (and is
// thus like a std::map or std::unordered_map. See
// https://stackoverflow.com/a/35293958/4081973 for more details.
template< typename T, typename = void > struct IsMappish : std::false_type { };

template<typename T> struct IsMappish< T, std::void_t< typename T::key_type,
    typename T::mapped_type, decltype(
    std::declval< T& >()[ std::declval< const typename T::key_type& >() ] )
    > > : std::true_type { };

/// Convert a JSON object to multiple data types using a common interface
template < typename T > bool convert_json( const marley::JSON& json,
  T& result )
{
  // Allow a trivial conversion to JSON itself (handy for populating a
  // container of JSON objects via recursive use of this function)
  if constexpr ( std::is_same_v< marley::JSON, T > ) {
    result = json;
    return true;
  }

  // This option includes both integers and boolean types, so we need to
  // explicitly distinguish between them
  if constexpr ( std::is_integral_v< T > ) {

    // If a boolean was requested just retrieve it and return
    if constexpr ( std::is_same_v< bool, T > ) {
      result = json.to_bool_or_throw();
      return true;
    }

    // For other integer types, first retrieve a value as a long int (the
    // standard representation for marley::JSON integers)
    long temp_long = json.to_long_or_throw();

    // If the requested type is unsigned, then double-check that we're not
    // working with a negative value stored in the JSON object. If we are,
    // then complain
    if constexpr ( std::is_unsigned_v< T > ) {
      if ( temp_long < 0 ) throw marley::Error( "Negative value "
        + std::to_string(temp_long) + " encountered when retrieving an"
        " unsigned integer from the JSON object '" + json.dump_string() + "'" );
    }

    // Things seem ok, so cast to the desired output type before returning
    result = static_cast< T >( temp_long );
    return true;
  }

  // For floating-point types, first retrieve a value as a double (the
  // standard representation for marley::JSON floating-point numbers)
  if constexpr ( std::is_floating_point_v< T > ) {
    double temp_double = json.to_double_or_throw();

    // Do any needed type conversion (e.g., to float) via this cast, then
    // return
    result = static_cast< T >( temp_double );
    return true;
  }

  if constexpr ( std::is_same_v< std::string, T > ) {
    result = json.to_string_or_throw();
    return true;
  }

  // If we've been handed a container that implements the push_back() and
  // clear() methods, then we will assume a JSON array is meant to be parsed.
  // Note that although JSON arrays are allowed to contain elements of
  // multiple data types, this implementation assumes that they are all
  // representable as T::value_type.
  if constexpr ( HasPushBack< T >::value && HasClear< T >::value ) {

    // Empty the container's existing contents
    result.clear();

    // Complain if we're not actually working with a JSON array. Default
    // to using the input JSON object itself if we weren't handed a key.
    if ( !json.is_array() ) throw marley::Error( "Invalid JSON array '"
      + json.dump_string() + "'" );

    // Fill the container with the JSON array elements using a dummy object
    // to allow for recursion
    auto elements = json.array_range();
    for ( const auto& el : elements ) {

      typename T::value_type temp_val;
      bool ok = convert_json( el, temp_val );

      // Save the element in the container if parsing went all right
      if ( ok ) result.push_back( temp_val );
      // Otherwise, complain
      else throw marley::Error( "Invalid array entry '" + el.dump_string()
        + "' found when parsing JSON array '" + json.dump_string() + "'" );
    }

    return true;
  }

  // If we've been handed a container that is "mappish" and implements the
  // clear() method, then we will assume a JSON object is meant to be parsed
  // into a map of key-value pairs. Note that this implementation requires
  // elements of the JSON object to have the same data type.
  if constexpr ( IsMappish< T >::value && HasClear< T >::value ) {

    // Empty the container's existing contents
    result.clear();

    // Complain if we're not actually working with a JSON object. Default
    // to using the input JSON object itself if we weren't handed a key.
    if ( !json.is_object() ) throw marley::Error( "Invalid JSON object '"
      + json.dump_string() + "'" );

    // Fill the container with the JSON object elements using recursion
    auto elements = json.object_range();
    for ( const auto& el : elements ) {

      std::string element_key = el.first;
      const marley::JSON& element_json = el.second;

      typename T::mapped_type temp_val;
      bool ok = convert_json( element_json, temp_val );

      // Save the element in the container if parsing went all right
      if ( ok ) result[ element_key ] = temp_val;
      // Otherwise, complain
      else throw marley::Error( "Invalid map entry '"
        + element_json.dump_string() + "' found when parsing JSON object '"
        + json.dump_string() + "'" );
    }

    return true;
  }

  // If we get here, then no conversion is implemented for the requested data
  // type
  return false;
}

// Alternate version that provides a default value if conversion fails
template < typename T > bool convert_json( const marley::JSON& json,
  T& result, T def_val )
{
  bool ok = convert_json( json, result );
  if ( !ok ) result = def_val;
  return ok;
}

// Alternate version that returns the converted value or a default one while
// storing the boolean flag indicating success (true) or failure (false)
template < typename T > T assign_from_json( const marley::JSON& json,
  bool& ok, T def_val = T() )
{
  T temp;
  ok = convert_json( json, temp, def_val );
  return temp;
}

/// Interface that converts an element of a JSON object with a given key
template < typename T > bool get_from_json( const std::string& key,
  const marley::JSON& json, T& result )
{
  json.check_if_object( key );
  if ( !json.has_key(key) ) return false;
  const marley::JSON& element = json.at( key );
  return convert_json( element, result );
}

// Alternate version that provides a default value if retrieval fails
template < typename T > bool get_from_json( const std::string& key,
  const marley::JSON& json, T& result, T def_val )
{
  bool ok = get_from_json( key, json, result );
  if ( !ok ) result = def_val;
  return ok;
}

// Alternate version that returns the retrieved value or a default one while
// storing the boolean flag indicating success (true) or failure (false)
template < typename T > T assign_from_json( const std::string& key,
  const marley::JSON& json, bool& ok, T def_val = T() )
{
  T temp;
  ok = get_from_json( key, json, temp, def_val );
  return temp;
}
