// -------------------------------------------------------------- -*- C++ -*-
// File: ilomipex1.cpp
// Version 12.6.3
// --------------------------------------------------------------------------
// Licensed Materials - Property of IBM
// 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
// Copyright IBM Corporation 2000, 2015. All Rights Reserved.
//
// US Government Users Restricted Rights - Use, duplication or
// disclosure restricted by GSA ADP Schedule Contract with
// IBM Corp.
// --------------------------------------------------------------------------
//
// ilomipex1.cpp - Entering and optimizing a MIP problem

#include <ilcplex/ilocplex.h>
#include <string.h>
#include <sstream>
#include <algorithm>
#include <unordered_map>

ILOSTLBEGIN

#define EPSZERO        1.0E-10

static void
   populatebyrow(IloModel model, IloNumVarArray var, IloRangeArray con);

int find(string arr[], int len, string seek)
{
  for (int i = 0; i < len; ++i)
    {
      if (arr[i] == seek) return i;
    }
  cout << "not found" << seek << endl;
  return -1;
}

int
main (int argc, char **argv) {
   IloEnv env;
   try {
      IloModel model(env);
      IloCplex cplex(env);
      // IloCplex::importModel(model,"./maxlength.lp")
      IloObjective   obj;
      IloNumVarArray var(env);
      IloRangeArray rng(env);
      cplex.importModel(model, argv[1], obj, var, rng);
      cplex.extract(model);
      
      cplex.setParam(IloCplex::EpRHS, 1e-9);
      cplex.setParam(IloCplex::EpOpt, 1e-9);
      cplex.setParam(IloCplex::EpInt, 0);
      cplex.setParam(IloCplex::VarSel, 3);
      cplex.setParam(IloCplex::ScaInd, 1);
      cplex.setParam(IloCplex::BndStrenInd, 1);
      cplex.setParam(IloCplex::NodeFileInd, 2);      


      cplex.solve();
      
      env.out() << "Solution status = " << cplex.getStatus() << endl;
      env.out() << "Solution value  = " << cplex.getObjValue() << endl;

      IloNumArray vals(env);
      cplex.getValues(vals, var);

      int count = 0;
      int iterN = 400;
      vector<string> solns {};
      while(count < iterN) {
      	int xIndex;
      	int yIndex;

      	for(int m=0; m< 6081; m=m+1) {
      	  if (vals[m] > 0.9 ){
	    string xname(var[m].getName());
	    bool found = (std::find(solns.begin(), solns.end(), xname) != solns.end());
	    if (!found) {
	      solns.push_back(xname);
	      xIndex=m;
	      printf("start: %d\n",xIndex);
	      env.out() << var[xIndex].getName() << endl;
	    }
	  }
      	}
      	for(int n=6081; n<12162; n=n+1) {
      	  if (vals[n] > 0.9){
	    string yname(var[n].getName());
	    bool found = (std::find(solns.begin(), solns.end(), yname) != solns.end());
	    if (!found) {
	      solns.push_back(yname);
	      yIndex=n;
	      printf("end: %d \n",yIndex - 4341);
	      env.out() << var[yIndex].getName() << endl;
	    }
      	  }
      	}
	env.out() << "-----------------------------------------------------" << endl;
	for (auto i: solns)
	  env.out() << i << ' ';
	env.out() << endl;
	env.out() << "-----------------------------------------------------" << endl;

	for(int i = 0; i < cplex.getNrows() - count*2; i=i+1) {
	  string name(rng[i].getName());
	  if (name == "start" || name == "end" || name == "allow_one_insertion") {
	    IloNum LB = (IloNum)count + 2.0;
	    IloNum UP = (IloNum)count + 2.0;
	    rng[i].setBounds(LB,UP);
	    env.out() << rng[i].getLB() << " and " << rng[i].getUB() << endl;
	  }
	}

	IloRangeArray c(env);
	c.add(var[xIndex] == 1);
	c.add(var[yIndex] == 1);
	c[0].setName("solx");
	c[1].setName("soly");
	model.add(c);

      	cplex.solve();
      	count = count + 1;
      	env.out() << "Solution status = " << cplex.getStatus() << endl;
      	env.out() << "Solution value  = " << cplex.getObjValue() << endl;
	cplex.getValues(vals, var);
      }
   }
   catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }

   env.end();
   return 0;

}  // END main
