#pragma once

struct TimeMeasure {
	std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
};

inline TimeMeasure StartTimeMeasure(){
	TimeMeasure Time = {};
	Time.start_time = std::chrono::high_resolution_clock::now();
	return Time;
}

inline void LogElapsedTime(const TimeMeasure& Time){
	auto finish_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_time = finish_time - Time.start_time;
	LogSystem::write("");
	LogSystem::write("->-> Total Elapsed Time: " + std::to_string(elapsed_time.count()) + " s");
	LogSystem::write("");
}
