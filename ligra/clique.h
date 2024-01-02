#ifndef CLIQUE_H
#define CLIQUE_H

#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include "parallel.h"
#include "quickSort.h"
#include "utils.h"


using namespace std;



struct edge_list{
	uintE n;
	uintE* edges; // out-edges and out-degrees and their offset  
	uintT* offsets;
	uintT* Degrees;

	bool Asymmetric; //is directed?
	uintE* inEdges; // in-edges and in-degrees and their offset
	uintT* Toffsets;
	uintT* DegreesT;

	void del() {
		free(edges);
		free(offsets);
		free(Degrees);
		if (Asymmetric) {
			free(inEdges);
			free(Toffsets);
			free(DegreesT);
		}
	}

	edge_list(int N, uintE* E, uintT* O, uintE* inE, uintT* inO) {
		n = N;
		edges = E;
		offsets = O;
		Degrees = newA(uintT,n);
		DegreesT = newA(uintT,n);
		Asymmetric = false;
		{parallel_for (long i = 0; i < n; i ++) {
			uintT o = offsets[i];
      		intT d = offsets[i+1] - o;
    		if(d < 0 || d > n) { 
      			cout << "degree out of bounds: vertex "<<i<< " has out-degree "<< d<<endl; 
      			abort(); }
      		Degrees[i]=DegreesT[i] = d;
		}}
		if (inE != NULL && inO != NULL) {
			Asymmetric = true;
			inEdges = inO;
			Toffsets = inO;
			{parallel_for (long i = 0; i < n; i ++) {
				uintT o = Toffsets[i];
      			intT d = Toffsets[i+1] - o;
    			if(d < 0 || d > n) { 
      				cout << "degree out of bounds: vertex "<<i<< " has out-degree "<< d<<endl; 
      				abort(); }
      			DegreesT[i] = d;
			}}
			}
	}
	

	void sort() {
		{parallel_for (long i=0; i < n; i++) {
      		uintT o = offsets[i];
      		intT d = Degrees[i];
    		if(d > 0) {
      			quickSort(edges+o, d, less<intE>());
    		}
   		}}
   		if (Asymmetric) {
   			{parallel_for (long i=0; i < n; i++) {
      			uintT o = Toffsets[i];
      			intT d = DegreesT[i+1] - o;
    			if(d > 0) {
      				quickSort(inEdges+o, d, less<intE>());
    		}
   		}}
   		}
	}

	void encode (uintE* mask, uintE* Map, uintE* MapT) {
		if (!Asymmetric) encodeSym(mask,Map,MapT);
	}

	void encodeSym(uintE* mask, uintE* Map, uintE* MapT) {
		uintE nu = 0, sn = 0;
		uintE* NewEdges = newA(uintE, offsets[n]);
		uintE* NewOffsets = newA(uintE, n+1);
		uintE* NewDegrees = newA(uintE, n);
		uintE index = 0;
		for(long i = 0; i < n; i++) {
			uintT vi = Map[i];
			NewOffsets[i] = index;
			uintT o = offsets[vi];
			if (mask[vi] == -1) {
				for (long j = 0; j < Degrees[vi]; j ++) {
					NewEdges[index++] = MapT[edges[o+j]];
					nu++;
				}
			} else {
				for (long j = 0; j < Degrees[vi]; j ++) {
					if (mask[vi] != mask[edges[o+j]])
						NewEdges[index++] = MapT[edges[o+j]];
					else sn++;
					nu ++;
				}
				NewEdges[index++] = mask[vi] + n;
				nu --;
			}
			NewDegrees[i] = index - NewOffsets[i];
		}
		NewOffsets[n] = index;
		edges  = NewEdges;
		offsets = NewOffsets;
		Degrees = NewDegrees;
		cout<<"compressing complete. Saving edges:" << sn<< "total edge num :"<<nu<<endl;
	}

	void encodeAsym(uintE* mask) {
		uintE* NewEdges = newA(uintE, offsets[n]);
		uintE* NewOffsets = newA(uintE, n+1);
		uintE* NewDegrees = newA(uintE, n);
		uintE* NewInEdges = newA(uintE, offsets[n]);
		uintE* NewInOffsets = newA(uintE, n+1);
		uintE* NewDegreesT = newA(uintE, n);
		uintE index = 0;
	}
};	



using Clique_Node = pair<uintT, vector<uintE>>;

struct clique{
	set<Clique_Node> C;
	uintE* Clique_index;
	edge_list* graph;
	uintE* Map,* MapT;
	uintE* mask;
	long max = 3;

	void del() {
		free(Clique_index);
		free(Map);
		free(MapT);
		free(mask);
		graph->del();
	}

	void encode(int N, uintE* E, uintT* O) {
		graph = new edge_list(N, E, O, NULL, NULL);
		Map = newA(uintE,N);
		MapT = newA(uintE,N);
		mask = newA(uintE,N);
		for (int i = 0; i < N; i++) mask[i] = -1;
		find_max_clique();
		ReMap();
		graph->encode(mask,Map,MapT);
	}

	void encode(int N, uintE* E, uintT* O, uintE* inE, uintT* inO) {
		graph = new edge_list(N, E, O, inE, inO);
	}

	bool MacCliqueHeu(Clique_Node& result, vector<uintE> S1) {
		if (S1.size() == 0) {
			if (result.second.size() > max) return true;
			else return false;
		}
		vector<uintE> S2, S3;
		long Num = S1.size();
		uintT MaxOutIndex = 0;
		for (int i = 0; i < Num; i ++) {
			if (graph->Degrees[S1[i]] > graph->Degrees[S1[MaxOutIndex]]) MaxOutIndex = i;
		}
		uintT MNode =  S1[MaxOutIndex];
		result.second.push_back(MNode);
		mask[MNode] = result.first;
		S1.erase(S1.begin() + MaxOutIndex);
		uintE o = graph->offsets[MNode];
		for (int i = 0; i < graph->Degrees[MNode]; i ++) {
			uint node = graph->edges[o+i];
			if (mask[node] == -1)
				S2.push_back(node);
		} 
		set_intersection(S1.begin(),S1.end(),S2.begin(),S2.end(),back_inserter(S3));
		return MacCliqueHeu(result,S3);
	}

	uintT cliqueHEU(uintT j, Clique_Node& result){
		vector<uintE> S1;
		uintT i;
		for (i = j; i < graph->n; i ++) {
			if (mask[i]== -1 && graph->Degrees[i] >= max) {
				mask[i] = result.first;
				result.second.push_back(i);
				for (int m = 0; m < graph->Degrees[i]; m ++) {
					uintT o = graph->offsets[i] + m;
					uintE nrr = graph->edges[o];
					if (mask[nrr] == -1){
						S1.push_back(nrr);
					}
				}
				if(MacCliqueHeu(result, S1)) break;
				else {
					for (int k = 0; k < result.second.size(); k ++) {
						mask[result.second[k]] = -1;
					}
					result.second.clear();
					S1.clear();
				}
				}
		}
		return i;
	}

	void find_max_clique() {
		uintT num = 0;
		uintT lastFind = 0;
		while(true) {
			Clique_Node result;
			result.first = num++;
			lastFind = cliqueHEU(lastFind, result);

			if (result.second.size()>=3) C.insert(result);
			if (lastFind >= graph->n) break;
		}
	}

	void ReMap(){
		Map = newA(uintE,graph->n);
		MapT = newA(uintE,graph->n);
		Clique_index = newA(uintE,C.size());
		uintE u = graph->n - 1;
		for (set<Clique_Node>::iterator it = C.begin(); it != C.end(); it++) {
			for (int j = 0; j < it->second.size(); j ++) {
				cout<<"Map "<<u<<" -:"<<it->second[j] << " ";
				Map[u--] = it->second[j];
			}
			Clique_index[it->first] = u+1;
			cout<<endl;
		}

		for (int i = graph->n - 1; i >= 0 ; i --) {
			if (mask[i] == -1) {
				Map[u--] = i;
			}
		}
		for (int i = 0; i < graph->n; i++) {
			MapT[Map[i]] = i;
		}
		// view_Arr(Map,graph->n);
	}

	void view_clique() {
		cout<<"Clique number :"<<C.size()<<endl;
		uintE nu=0;
		for (set<Clique_Node>::iterator it = C.begin(); it != C.end(); it++) {
			// cout<<"clique "<<it->first<<"  ";
			for (int j = 0; j < it->second.size(); j ++) {
				// cout<<it->second[j]<<" ";
				nu ++;
			}
			// cout<<endl;
		}
		cout<<"averg:"<< nu*1.0 / C.size()<<endl;
	}

};

#endif