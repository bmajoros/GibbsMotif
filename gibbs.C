/****************************************************************
 gibbs.C
 william.majoros@duke.edu

 This is open-source software, governed by the ARTISTIC LICENSE 
 (see www.opensource.org).
 ****************************************************************/
#include <math.h>
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/Sequence.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Constants.H"
#include "Motif.H"
using namespace std;
using namespace BOOM;

const float PSEUDOCOUNT=0.1;

enum SamplerType {
  GIBBS,
  METROPOLIS_HASTINGS
};

class Application {
  int L; // motif length
  PureDnaAlphabet alphabet;
  Vector<Sequence> seqs;
  SamplerType samplerType;
  float MH_maxProposalSize;
  void load(const String &filename);
  int samplePos(int seqIndex,Array1D<int> &positions);
  Motif *trainMotif(Array1D<int> &positions,int omit=-1);
  int sampleGibbs(int seqIndex,Array1D<int> &positions);
  int sampleMH(int seqIndex,Array1D<int> &positions);
  void getPositionDistr(Array1D<float> &P,Motif &,Sequence &);
  int samplePosition(Array1D<float> &P);
  void downSampleArray(const Array1D<float> &P,Array1D<float> &Q,
		       Array1D<int> &q);
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
  CommandLine cmd(argc,argv,"Hp:");
  if(cmd.numArgs()!=4)
    throw String("\n\
gibbs [options] <*.fasta> <length> <iterations> <#restarts>\n\
      -H = use Metropolis-Hastings instead of Gibbs\n\
      -p N = MH proposal considers at most N positions when choosing\n\
             (default N=1 implies proposal is uniform random function,\n\
              N=0 implies no limit => proposal is fully informed)\n\
");
  const String fasta=cmd.arg(0);
  L=cmd.arg(1).asInt();
  const int iter=cmd.arg(2).asInt();
  const int restarts=cmd.arg(3).asInt();
  samplerType=cmd.option('H') ? METROPOLIS_HASTINGS : GIBBS;
  MH_maxProposalSize=cmd.option('p') ? cmd.optParm('p').asInt() : 1;
  if(MH_maxProposalSize<1) MH_maxProposalSize=INT_MAX;

  // Load sequences
  load(fasta);
  int N=seqs.size(); 
  if(N<1) throw "no sequences loaded";

  // Perform sampling iterations
  Motif *bestM=NULL;
  float bestLL=NEGATIVE_INFINITY;
  for(int k=0 ; k<restarts ; ++k) {
    Array1D<int> Z(N);
    for(int i=0 ; i<N ; ++i) Z[i]=RandomNumber(0,seqs[i].getLength()-L);
    for(int i=0 ; i<iter ; ++i)
      for(int j=0 ; j<N ; ++j)
	Z[j]=samplePos(j,Z);
    Motif *M=trainMotif(Z);
    float LL=M->likelihood(seqs,Z);
    if(LL>bestLL) { bestM=M; bestLL=LL; }
    else delete M;
  }

  // Compile motif & emit
  cout<<*bestM<<endl;
  bestM->getConsensus().printOn(cout,alphabet);cout<<endl;
  delete bestM;

  return 0;
}



void Application::load(const String &filename)
{
  FastaReader reader(filename,alphabet);
  String def, seq;
  while(reader.nextSequence(def,seq)) {
    seqs.push_back(Sequence(seq,alphabet));
  }
}



Motif *Application::trainMotif(Array1D<int> &positions,int omit)
{
  Motif &M=*new Motif(L);
  M.setAllTo(PSEUDOCOUNT);
  int N=seqs.size();
  for(int i=0 ; i<N ; ++i) {
    if(i==omit) continue;
    Sequence &seq=seqs[i];
    const int begin=positions[i];
    for(int j=begin ; j<begin+L ; ++j)
      ++M.index(j-begin,seq[j]);
  }
  M.normalize();
  return &M;
}



int Application::samplePos(int seqIndex,Array1D<int> &positions)
{
  switch(samplerType)
    {
    case GIBBS: return sampleGibbs(seqIndex,positions);
    case METROPOLIS_HASTINGS: return sampleMH(seqIndex,positions);
    }
  INTERNAL_ERROR;
}



int Application::sampleGibbs(int seqIndex,Array1D<int> &positions)
{
  // Build motif model
  Motif &M=*trainMotif(positions,seqIndex);

  // Compute sampling distribution for positions in target sequence
  Sequence &seq=seqs[seqIndex];
  int len=seq.getLength();
  int end=len-L;
  Array1D<float> P(end+1);
  getPositionDistr(P,M,seq);

  // Sample position
  int pos=samplePosition(P);

  // clean up
  delete &M;
  return pos;
}



int Application::samplePosition(Array1D<float> &P)
{
  float r=Random0to1();
  float sum=0, pos;
  int end=P.size();
  for(pos=0 ; pos<end ; ++pos) {
    sum+=P[pos];
    if(r<=sum) break;
  }
  return pos;
}



void Application::getPositionDistr(Array1D<float> &P,Motif &M,Sequence &seq)
{
  const int len=P.size();
  for(int i=0 ; i<len ; ++i) P[i]=exp(M.likelihood(seq,i));
  float sum=0;
  for(int i=0 ; i<len ; ++i) sum+=P[i];
  for(int i=0 ; i<len ; ++i) P[i]/=sum;
}



void Application::downSampleArray(const Array1D<float> &P,Array1D<float> &Q,
				  Array1D<int> &q)
{
  int n=P.size(), m=Q.size();
  Array1D<bool> hit(n); hit.setAllTo(false);
  for(int i=0 ; i<m ; ++i) {
    int pos=RandomNumber(n);
    while(hit[pos]) pos=(pos+1)%n;
    Q[i]=P[pos]; q[i]=pos;
    hit[pos]=true;
  }
  float sum=0;
  for(int i=0 ; i<m ; ++i) sum+=Q[i];
  for(int i=0 ; i<m ; ++i) Q[i]/=sum;
}



int Application::sampleMH(int seqIndex,Array1D<int> &positions)
{
  Sequence &seq=seqs[seqIndex];
  const int len=seq.getLength();
  Motif &M=*trainMotif(positions,seqIndex);
  const int prevPos=positions[seqIndex];
  const L=M.length();
  int last=len-L;
  float prior=log(1.0/(last+1));
  float lik1=M.likelihood(seqs,positions);

  // Until acceptance...
  int newPos;
  int iter=1;
  float oldProposal=1.0/(last+1), newProposal=1.0/(last+1);
  while(1) {
    // Make a proposal
    Array1D<float> P(last+1);
    getPositionDistr(P,M,seq);
    const int k=MH_maxProposalSize<P.size() ? MH_maxProposalSize : P.size();
    Array1D<float> Q(k); Array1D<int> q(k);
    downSampleArray(P,Q,q);
    int j=samplePosition(Q);
    newPos=q[j];
    //cout<<P<<endl<<Q<<endl<<"j="<<j<<" newPos="<<newPos<<endl;
    oldProposal=P[prevPos]; newProposal=P[newPos];

    // Evaluate Hastings ratio
    positions[seqIndex]=newPos;
    Motif &m=*trainMotif(positions);
    float lik2=m.likelihood(seqs,positions);
    float prior1=prior, prior2=prior;
    float pi1=lik1+prior1, pi2=lik2+prior2;
    float hastings=exp(pi2-pi1+log(oldProposal)-log(newProposal));
    //cout<<"H="<<hastings<<" lik1="<<lik1<<" lik2="<<lik2<<" prior1="<<prior1<<" prior2="<<prior2<<endl;
    if(!isFinite(hastings)) {
      cout<<M<<endl; cout<<m<<endl;
      cout<<positions<<endl;
      cout<<"L="<<L<<endl;
      throw String("hastings=")+hastings+" lik1="+lik1+" lik2="+lik2+" prior1="+prior1+" prior2="+prior2+" len="+len+" prevPos="+prevPos+" newPos="+newPos+" seqIndex="+seqIndex;
      delete &m;
    }
    // Hastings Ratio = pi(x2)p(x1|x2) / pi(x1)p(x2|x1)

    // Accept or reject
    float accept=fmin(1.0,hastings);
    float r=Random0to1();
    if(r<=accept) break;
    if(++iter>1000) {
      cerr<<"WARNING!  too many iterations, breaking..."<<endl;
      break; }
  }
  delete &M;
  positions[seqIndex]=prevPos;
  return newPos;
}



