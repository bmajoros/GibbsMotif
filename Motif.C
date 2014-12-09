/****************************************************************
 Motif.C
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <math.h>
#include <iostream>
#include <fstream>
#include "BOOM/Random.H"
#include "Motif.H"
using namespace std;
using namespace BOOM;

Motif::Motif(int length)
{
  if(length>0) pwm.resize(length,alphabet.size());
}



Motif::Motif(const String &filename)
{
  const int nAlpha=alphabet.size();
  ifstream is(filename.c_str());
  if(!is.good()) throw String("Cannot open file ")+filename;
  int L;
  is>>L;
  resize(L);
  for(int i=0 ; i<L ; ++i) {
    for(int j=0 ; j<nAlpha ; ++j) {
      float f;
      is>>f;
      pwm[i][j]=f;
    }
  }
}



Alphabet &Motif::getAlphabet()
{
  return alphabet;
}



void Motif::save(const String &filename)
{
  ofstream os(filename.c_str());
  int L=length(), N=alphabet.size();
  os<<L<<endl;
  for(int i=0 ; i<L ; ++i) {
    for(int j=0 ; j<N ; ++j) os<<pwm[i][j]<<"\t";
    os<<endl;
  }
}



int Motif::length()
{
  return pwm.getFirstDim();
}



void Motif::resize(int length)
{
  pwm.resize(length,alphabet.size());
}



float &Motif::index(int pos,Symbol s)
{
  return pwm[pos][s];
}



Symbol Motif::sample(int pos)
{
  float r=Random0to1(), sum=0;
  int n=alphabet.size();
  for(int i=0 ; i<n ; ++i) {
    sum+=pwm[pos][i];
    if(r<=sum) return i;
  }
  return n-1;
}



void Motif::addPseudocount(float c)
{
  int L=length(), n=alphabet.size();
  for(int i=0 ; i<L ; ++i)
    for(int j=0 ; j<n ; ++j)
      pwm[i][j]+=c;
}



void Motif::normalize()
{
  int L=length(), n=alphabet.size();
  for(int i=0 ; i<L ; ++i) {
    float sum=0;
    for(int j=0 ; j<n ; ++j) sum+=pwm[i][j];
    for(int j=0 ; j<n ; ++j) pwm[i][j]/=sum;
  }
}



void Motif::printOn(ostream &os) const
{
  int L=length(), N=alphabet.size();
  os<<L<<endl;
  for(int i=0 ; i<L ; ++i) {
    for(int j=0 ; j<N ; ++j) os<<pwm[i][j]<<"\t";
    os<<endl;
  }
}



ostream &operator<<(ostream &os,const Motif &M)
{
  M.printOn(os);
  return os;
}



float Motif::likelihood(Sequence &S,int pos)
{
  int L=length();
  float P=0;
  for(int i=0 ; i<L ; ++i)
    P+=log(pwm[i][S[pos+i]]);
  return P;
}



float Motif::likelihood(Vector<Sequence> &A,Array1D<int> &positions)
{
  int n=A.size();
  float sum=0;
  for(int i=0 ; i<n ; ++i) sum+=likelihood(A[i],positions[i]);
  return sum;
}



void Motif::setAllTo(float f)
{
  pwm.setAllTo(f);
}



Sequence Motif::getConsensus()
{
  const int L=length(), n=alphabet.size();
  Sequence S;
  for(int i=0 ; i<L ; ++i) {
    int max=0;
    for(int j=0 ; j<n ; ++j) if(pwm[i][j]>pwm[i][max]) max=j;
    S+=Symbol(max);
  }
  return S;
}



