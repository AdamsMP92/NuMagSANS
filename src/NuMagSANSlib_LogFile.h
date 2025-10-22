// File         : NuMagSANSlib_LogFile.h
// Author       : Dr. Michael Philipp ADAMS
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 16 October 2025
// OS           : Linux Ubuntu
// Language     : CUDA C++

#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>

namespace LogSystem {

inline std::ofstream logFile;

// open log file
inline void initLog(const std::string& filename = "log.txt")
{
    logFile.open(filename, std::ios::out | std::ios::app);
    if (!logFile.is_open()) {
        std::cerr << "Could not open log file: " << filename << std::endl;
        return;
    }

    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    logFile << "# Log started: " << std::put_time(std::localtime(&now), "%F %T") << "\n";
    logFile.flush();
    std::cout << "Logging to: " << filename << std::endl;
}

// write message to log file
inline void write(const std::string& msg)
{
    if (!logFile.is_open()) return;

    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    logFile << std::put_time(std::localtime(&now), "%F %T") << " | " << msg << '\n';
    logFile.flush();

    // optional: auch aufs Terminal
    std::cout << std::put_time(std::localtime(&now), "%F %T") << " | " << msg << std::endl;
}

// close log file
inline void close()
{
	auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    if (logFile.is_open()) {
        logFile << "# Log closed:" << std::put_time(std::localtime(&now), "%F %T") << "\n";
        logFile.close();
    }
}

} // namespace LogSystem
