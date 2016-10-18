/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMacIonize is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file Log.hpp
 *
 * @brief General interface for logging classes.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef LOG_HPP
#define LOG_HPP

#include <ctime>
#include <sstream>
#include <string>

/**
 * @brief Possible logging levels.
 *
 * The higher the level of a log message, the more likely it is to be written to
 * the log.
 */
enum LogLevel {
  LOGLEVEL_INFO = 0,
  LOGLEVEL_STATUS,
  LOGLEVEL_WARNING,
  LOGLEVEL_ERROR
};

/**
 * @brief General interface for loggers.
 *
 * This abstract class implements log routines for different log levels, but
 * does not implement the actual writing of log messages, which needs to be done
 * in a child class.
 */
class Log {
private:
  /*! @brief Log level. */
  LogLevel _level;

  /**
   * @brief Get a time stamp to append to log messages.
   *
   * The time stamp has the format HH:MM:SS.
   *
   * @return Time stamp.
   */
  inline std::string get_timestamp() {
    std::time_t timestamp = std::time(NULL);
    std::tm *time = std::localtime(&timestamp);
    std::stringstream timestream;
    timestream << time->tm_hour << ":" << time->tm_min << ":" << time->tm_sec;
    return timestream.str();
  }

protected:
  /**
   * @brief Write the given message to the log.
   *
   * This method should be implemented in a child class.
   *
   * @param message Message to write to the log.
   */
  virtual void write_message(std::string message) = 0;

  /**
   * @brief Get a description for the given LogLevel that can be written to the
   * logger.
   *
   * We make this method virtual, so that we can use colours in the terminal
   * implementation.
   *
   * @param level LogLevel of the message.
   */
  virtual std::string get_levelname(LogLevel level) {
    switch (level) {
    case LOGLEVEL_INFO:
      return "info";
    case LOGLEVEL_STATUS:
      return "status";
    case LOGLEVEL_WARNING:
      return "warning";
    case LOGLEVEL_ERROR:
      return "error";
    default:
      return "message";
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param level LogLevel of the logger.
   */
  inline Log(LogLevel level) : _level(level) {}

  /**
   * @brief Write the given message with the given log level to the logger.
   *
   * The message is only written out if the given log level is higher or equal
   * to the logger level.
   *
   * @param level LogLevel of the message.
   * @param message Message.
   */
  inline void write_message(LogLevel level, std::string message) {
    if (level >= _level) {
      std::stringstream messagestream;
      messagestream << get_timestamp() << ": " << get_levelname(level) << ": "
                    << message;
      write_message(messagestream.str());
    }
  }

  /**
   * @brief Convert the given template argument to a std::string by piping it
   * through a std::stringstream.
   *
   * @param t Template argument.
   * @return std::string containing a representation of the argument.
   */
  template < typename T > std::string get_message(T t) {
    std::stringstream stream;
    stream << t;
    return stream.str();
  }

  /**
   * @brief Recursive variadic version of get_message() that takes a variable
   * number of arguments and parses them one by one.
   *
   * @param t Next argument in the list.
   * @param args Remaining arguments.
   * @return std::string containing the arguments one after another.
   */
  template < typename T, typename... Args >
  std::string get_message(T t, Args... args) {
    std::stringstream stream;
    stream << t << get_message(args...);
    return stream.str();
  }

  /**
   * @brief Write an info message to the logger.
   *
   * Info messages have the lowest priority and are only written if the logger
   * has the lowest level.
   *
   * @param args Things to write to the logger.
   */
  template < typename... Args > inline void write_info(Args... args) {
    write_message(LOGLEVEL_INFO, get_message(args...));
  }

  /**
   * @brief Write a status message to the logger.
   *
   * Status messages are meant to reflect the current status of the program and
   * should be written out regularly, so that the user has an idea of what the
   * program is currently doing.
   *
   * @param args Things to write to the logger.
   */
  template < typename... Args > inline void write_status(Args... args) {
    write_message(LOGLEVEL_STATUS, get_message(args...));
  }

  /**
   * @brief Write a warning message to the logger.
   *
   * Warning messages indicate a problematic situation that can be overcome by
   * the program at run time, but should be noted by the user.
   *
   * @param args Things to write to the logger.
   */
  template < typename... Args > inline void write_warning(Args... args) {
    write_message(LOGLEVEL_WARNING, get_message(args...));
  }

  /**
   * @brief Write an error message to the logger.
   *
   * Error messages indicate a problematic situation that is fatal for the
   * program. These are usually handled by calling the error macro.
   *
   * @param args Things to write to the logger.
   */
  template < typename... Args > inline void write_error(Args... args) {
    write_message(LOGLEVEL_ERROR, get_message(args...));
  }
};

#endif // LOG_HPP
