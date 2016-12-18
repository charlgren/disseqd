#ifndef DISSEQD_CPP
#define DISSEQD_CPP

#define LOGSCALE

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <iostream>
#include "hmm.cpp"

using namespace std;

static void usage(string name)
{
  cerr<<"Usage: "<<name<<" -s seqfile (-m file / -i seqfile(s) -k kmer -n nstates) [-o prefix] [-r rounds]\n"
      <<"Options:\n"
      <<"\t-h,--help\tShow this message\n"
      <<"Mandatory:\n"
      <<"\t-s,--sequence seqfile\tFile with sequence(s) to parse\n"
      <<"Alternative:\n"
      <<"\t-m,--model file\t\tFile with hmm probabilities\n"
      <<"or\n"
      <<"\t-i,--in seqfile[,seqfiles]\tFile with sequence(s) for each state for training\n"
      <<"\t-k,--kmer int\t\tNumber of consequtive nucleotides to consider for training\n"
      <<"\t-n,--nstates int\tNumber of states in model for training\n"
      <<"Optional:\n"
      <<"\t-o,--out prefix\t\tPrefix for .model.txt and .decode.txt outfiles\n"
      <<"\t-r,--rounds int\t\tNumber of training iterations\n"
      <<endl;
}

int main(int argc, char **argv)
{
  //input parameters
  string infile = "";
  int kmer = 0;
  string modelfile = "";
  int nstates = 0;
  string outfile = "";
  int rounds = 0;
  string sequencefile = "";
  int c = 0;
//------------------------------- read arguments -------------------------------
  while(1)
  {
    static struct option long_options[]=
    {
      //these options set a flag
      {"debug", no_argument,  &debug_flag,  1},
      {"verbose", no_argument,  &verbose_flag,  1},
      //these options don't set a flag
      {"help", no_argument,  0,  'h'},
      {"in", required_argument,  0,  'i'},
      {"kmer", required_argument,  0,  'k'},
      {"model", required_argument,  0,  'm'},
      {"nstates", required_argument,  0,  'n'},
      {"out", required_argument,  0,  'o'},
      {"rounds", required_argument,  0,  'r'},
      {"sequence", required_argument,  0,  's'},
      {0, 0,  0,  0}
    };
    //getopt_long stores the option index here
    int option_index = 0;

    c = getopt_long(argc, argv, "hi:k:m:n:o:r:s:", long_options, &option_index);

    if(c==-1){break;}

    switch(c)
    {
      //when is case 0 entered?
      case 0:
        //if this option set a flag, do nothing else now
        if(long_options[option_index].flag != 0){break;}
        printf("option --%s", long_options[option_index].name);
        if(optarg){printf(" with arg %s", optarg);}
        printf("\n");
        break;
      case 'h':
        usage(argv[0]);
        exit(EXIT_FAILURE);
      case 'i':
        printf("option --in with value '%s'\n", optarg);
        infile = optarg;
        break;
      case 'k':
        printf("option --kmer with value '%s'\n", optarg);
        kmer = atoi(optarg);
        break;
      case 'm':
        printf("option --model with value '%s'\n", optarg);
        modelfile = optarg;
        break;
      case 'n':
        printf("option --nstates with value '%s'\n", optarg);
        nstates = atoi(optarg);
        break;
      case 'o':
        printf("option --out with value '%s'\n", optarg);
        outfile = optarg;
        break;
      case 'r':
        printf("option --rounds with value '%s'\n", optarg);
        rounds = atoi(optarg);
        break;
      case 's':
        printf("option --sequence with value '%s'\n", optarg);
        sequencefile = optarg;
        break;
      case '?':
        //getopt_long already printed an error message (invalid/unrecognized option)
        fprintf(stderr,"Usage: %s -s seqfile (-m modelfile / -i infile(s) -k kmer -n nstates) [-o outprefix] [-r rounds]\n", argv[0]);
        break;
      //when is default entered?
      default:
        abort();
    }
  }
//------------------------- handle remaining arguments -------------------------
  if(optind<argc)
  {
    printf("non-parsed non-option ARGV-elements:");
    while(optind<argc)
    {
      printf("\t%s", argv[optind++]);
    }
    putchar('\n');
  }
//------------------------------- run algorithms -------------------------------
  if(outfile=="") printf("missing --%s option, writing to stdout\n","outfile");
  if(sequencefile!="")
  {
    Hmm hmm = Hmm();
    string obs = hmm.read_fasta(sequencefile);
    //get model from file
    if(modelfile!="")
    {
      if(infile!="") printf("skipping --%s option\n", "in");
      if(kmer) printf("skipping --%s option\n", "kmer");
      if(nstates) printf("skipping --%s option\n", "nstates");
      hmm.read_model(modelfile);
    }
    //get model from parameters
    else if(infile!="" && kmer && nstates)
    {
      hmm = Hmm(nstates, kmer);
      vector<string> infiles = split(infile,',');
      if(infiles.size()!=nstates)
      {
        printf("please supply %d infiles\n",nstates);
      }
      for(int s=0; s<infiles.size(); s++)
      {
        hmm.read_emissions_fasta(s,infiles[s]);
      }
    }
    //missing arguments
    else
    {
      if(modelfile=="") printf("missing --%s option\nor alternatively:\n","model");
      if(infile=="") printf("missing --%s option\n", "in");
      if(!kmer) printf("missing --%s option\n", "kmer");
      if(!nstates) printf("missing --%s option\n", "nstates");
      exit(EXIT_FAILURE);
    }
    //train model
    if(rounds)
    {
      string done = "";
      for(int i=0; i<rounds; i++)
      {
        hmm.train(obs);
        done = hmm.decode(obs);
        if(done==""){rounds = i;}
        //hmm.decode(obs);
      }
      //write output
      hmm.write_model(outfile);
      hmm.write_decoding(outfile);
      printf("rounds: %d\n",rounds);
      exit(EXIT_SUCCESS);
    }
    //decode from model
    hmm.decode(obs);
    //write output
    hmm.write_model(outfile);
    hmm.write_decoding(outfile);
    //string str = "abc";
    //cout<<str.at(0)<<1./str.at(0)<<CHAR_BIT<<char(int(0.30*255))<<sizeof(char)<<str.size()<<probability(1.0/0)<<endl;
    exit(EXIT_SUCCESS);
  }
  else
  {
    printf("missing --%s option\n","sequence");
    exit(EXIT_FAILURE);
  }
}

#endif //DISSEQD_CPP
