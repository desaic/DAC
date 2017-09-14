#ifndef TIMER_HPP
#define TIMER_HPP
#include <chrono>

class Timer
{
public:
  Timer();
  void start();
  ///@return -1 if overflow happened
  int end();
  float getSeconds();
  float getMilliseconds();
  clock_t getClocks();

  void   startWall();
  void   endWall();
  double getSecondsWall();

private:
  clock_t t0, t1;

  std::chrono::time_point<std::chrono::system_clock> start_t;
  std::chrono::time_point<std::chrono::system_clock> finish_t;

};

#endif
