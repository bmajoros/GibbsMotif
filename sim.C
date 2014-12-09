/****************************************************************
 sim.C
 william.majoros@duke.edu

 This is open-source software, governed by the ARTISTIC LICENSE 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaWriter.H"
#include "Motif.H"
using namespace std;
using namespace BOOM;


class Application
{
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Application::Application()
{
  randomize();
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=5)
    throw String("sim <motif.txt> <#seqs> <length> <pseudocount> <fraction-present>");
  String motifFile=cmd.arg(0);
  int numSeqs=cmd.arg(1).asInt();
  int L=cmd.arg(2).asInt();
  float pseudocount=cmd.arg(3).asFloat();
  float presence=cmd.arg(4).asFloat();

  // Load motif
  Motif motif(motifFile);
  int motifLen=motif.length();
  Alphabet &alphabet=motif.getAlphabet();
  int nAlpha=alphabet.size();
  motif.addPseudocount(pseudocount);
  motif.normalize();

  // Generate sequences
  FastaWriter writer;
  for(int i=0 ; i<numSeqs ; ++i) {
    int pos=RandomNumber(0,L-motifLen);
    String seq(L,' ');
    for(int j=0 ; j<L ; ++j)
      seq[j]=alphabet.lookup(Symbol(RandomNumber(nAlpha)));
    if(Random0to1()<=presence)
      for(int j=pos ; j<pos+motifLen ; ++j)
	seq[j]=alphabet.lookup(motif.sample(j-pos));
    const int id=i+1;
    String def=String(">")+i+" /pos="+pos;
    writer.addToFasta(def,seq,cout);
  }

  return 0;
}

