#ifndef WINDOW_H
#define WINDOW_H

#include <deque>
#include <iostream>
#include <limits>
#include "signal.h"

using namespace std;

//piecewise-linear, right-continuous signals with area computation
class Window : public Signal {
 private:
  // pop the right most element
  void pop_back_update(void);
  // push an element from right
  void push_back_update(Sample, double);

 public:
  double area;
  
  Window() {}
  Window(Sample);

  double width();

// Move the left endpoint of the window to T.
  void updateFront(double, bool);
  void updateBackMax(Sample, bool);
  void updateBackMin(Sample, bool);
};

#endif
