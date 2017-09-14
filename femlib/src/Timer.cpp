#include "Timer.hpp"

#include<iostream>

Timer::Timer() :t0(0), t1(0)
{}

void Timer::start()
{
  startWall();
}

int Timer::end()
{
  endWall();
  return 0;
}

float Timer::getSeconds()
{
  std::chrono::duration<double> diff = finish_t - start_t;
  return diff.count();
}

float Timer::getMilliseconds()
{
  std::chrono::duration<double> diff = finish_t - start_t;
  return diff.count() / 1000;
}

clock_t Timer::getClocks()
{
  return t1 - t0;
}

void Timer::startWall()
{
  start_t = std::chrono::system_clock::now();
}

void Timer::endWall()
{
  finish_t = std::chrono::system_clock::now();
}

double Timer::getSecondsWall()
{
  std::chrono::duration<double> diff = finish_t - start_t;
  return diff.count();
}
