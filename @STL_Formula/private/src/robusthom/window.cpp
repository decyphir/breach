#include "window.h"

using namespace std;

Window::Window(Sample p){
  beginTime = p.time;
  endTime = p.time;
  area=0;
  push_front(p);
}

double Window::width(){
  return endTime - beginTime;
}

// pop the right most element and update the area
void Window::pop_back_update(){
  area -= back().area(endTime);
  endTime = back().time;
  pop_back();
  return;
}

// push an element from right and update the area
void Window::push_back_update(Sample newSample, double endT){
  push_back(newSample);
  endTime = endT;
  area += back().area(endTime);
  return;
}


// Move the left endpoint of the window to T.
void Window::updateFront(double T, bool end_candidate){
  Sample prevFront;
  double w = width();

  if(empty()){
    cout << "Illegal argument: Window must not be empty in updateFront" << endl;
    return;
  }
  
  // Dequeue old sample
  prevFront = front();
  pop_front();
  if(!end_candidate){
    push_front(Sample(T, prevFront.valueAt(T), prevFront.derivative));
  }

  beginTime = T;
  // Update area
  if(T > endTime){
    endTime = T+w;
    area = front().area(T+w);
  }else{
    area -= prevFront.area(T);
  }
}


// Fix the window from right to keep the value at t' to be max(t':endT).
void Window::updateBackMax(Sample newSample, bool new_candidate){
  bool cpFlag = false;
  double cpTime;
  Sample prevBack;

  if(empty()){
    beginTime = newSample.time;
    area = 0.;
    push_back_update(newSample, newSample.time);
    return;
  }
  // Pop all smaller elements
  while(!empty() && back().value <= newSample.value) {
    pop_back_update();
    cpFlag=true;
  }
  if(empty()){
    push_back_update(Sample(endTime, newSample.value, 0), newSample.time);
  }else if(cpFlag){
    if(back().derivative == 0){
      cpTime = endTime;
    }else{
      cpTime = fmin(endTime, back().timeIntersect(newSample.constant()));
    }
    prevBack = back();
    pop_back_update();
    push_back_update(prevBack, cpTime);
    push_back_update(Sample(cpTime, newSample.value, 0.), newSample.time);
  }else{
    prevBack = back();
    pop_back_update();
    push_back_update(prevBack, newSample.time);
  }

  if(new_candidate){
    push_back(newSample);
  }
}

// Fix the window from right to keep the value at t' to be min(t':endT).
void Window::updateBackMin(Sample newSample, bool new_candidate){
  bool cpFlag = false;
  double cpTime;
  Sample prevBack;

  if(empty()){
    beginTime = newSample.time;
    area = 0.;
    push_back_update(newSample, newSample.time);
    return;
  }
  // Pop all smaller elements
  while(!empty() && back().value >= newSample.value) {
    pop_back_update();
    cpFlag=true;
  }
  if(empty()){
    push_back_update(Sample(endTime, newSample.value, 0), newSample.time);
  }else if(cpFlag){
    if(back().derivative == 0){
      cpTime = endTime;
    }else{
      cpTime = fmin(endTime, back().timeIntersect(newSample.constant()));
    }
    prevBack = back();
    pop_back_update();
    push_back_update(prevBack, cpTime);
    push_back_update(Sample(cpTime, newSample.value, 0.), newSample.time);
  }else{
    prevBack = back();
    pop_back_update();
    push_back_update(prevBack, newSample.time);
  }

  if(new_candidate){
    push_back(newSample);
  }
}


