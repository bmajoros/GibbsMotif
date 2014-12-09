/****************************************************************
 Motif.H
 Copyright (C)2013 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Motif_H
#define INCL_Motif_H
#include <iostream>
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/Array2D.H"
#include "BOOM/Array1D.H"
#include "BOOM/Vector.H"
#include "BOOM/Sequence.H"
using namespace std;
using namespace BOOM;


class Motif {
  Array2D<float> pwm; // [pos,symbol]
  PureDnaAlphabet alphabet;
public:
  Motif(int length=0);
  void setAllTo(float);
  Motif(const String &filename);
  void save(const String &filename);
  int length();
  void resize(int length);
  float &index(int pos,Symbol);
  Alphabet &getAlphabet();
  Symbol sample(int pos);
  void addPseudocount(float);
  void normalize();
  void printOn(ostream &) const;
  float likelihood(Sequence &,int pos);
  float likelihood(Vector<Sequence> &,Array1D<int> &positions);
  Sequence getConsensus();
};

ostream &operator<<(ostream &,const Motif &);

#endif

