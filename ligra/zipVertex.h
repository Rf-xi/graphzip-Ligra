#ifndef ZIPVERTEX_H
#define ZIPVERTEX_H
#include "vertexSubset.h"
using namespace std;

namespace decode_zip {

  // Used by edgeMapDense. Callers ensure cond(v_id). For each vertex, decode
  // its in-edges, and check to see whether this neighbor is in the current
  // frontier, calling update if it is. If processing the edges sequentially,
  // break once !cond(v_id).
  template <class vertex, class F, class G, class VS>
  inline void decodeInNghBreakEarly(vertex* v, long v_id, VS& vertexSubset, F &f, G &g, bool parallel = 0) {
    long d = v->getInDegree();
    if (d==0) return;
    bool fl = false;
    if (!parallel || d < 1000) {
      for (size_t j=0; j<d-1; j++) {
        uintE ngh = v->getInNeighbor(j);
        if (vertexSubset.isIn(ngh)) {
          auto m = f.update(ngh, v_id);
          g(v_id, m);
        }
        if(!f.cond(v_id)) {fl=true;break;}
      }
      uintE ngh = v->getInNeighbor(d-1);
      if (!fl) {
        if (v->CliqueNode==NULL) {
          if (vertexSubset.isIn(ngh)) {
            auto m = f.update(ngh,v_id);
            g(v_id, m);
          }
        } else {
          if (ngh<v->CliqueNode->Tdegree) cout<<"error x."<<endl;
          for (long ki = v->CliqueNode->neighbors[0]; ki < v->CliqueNode->neighbors[1]; ki++) {
            if (ki == v_id) continue;
            if (vertexSubset.isIn(ki)) {
              auto m = f.update(ki, v_id);
              g(v_id,m);
            }
            if(!f.cond(v_id)) {break;}
          }
        }
      }
    } else {
      parallel_for(size_t j=0; j<d-1; j++) {
        uintE ngh = v->getInNeighbor(j);
        if (vertexSubset.isIn(ngh)) {
          auto m = f.updateAtomic(ngh, v_id);
          g(v_id, m);
        }
      }
      uintE ngh = v->getInNeighbor(d-1);
      if (v->CliqueNode==NULL) {
        if (vertexSubset.isIn(ngh)) {
          auto m = f.update(ngh,v_id);
          g(v_id, m);
        }
      } else { 
        if (ngh<v->CliqueNode->Tdegree) cout<<"error x."<<endl;
        parallel_for (long ki = v->CliqueNode->neighbors[0]; ki < v->CliqueNode->neighbors[1]; ki++) {
          if (ki!=v_id && vertexSubset.isIn(ki)) {
            auto m = f.update(ki, v_id);
            g(v_id, m);
          }
        }
      }
    }
  }

  // Used by edgeMapDenseForward. For each out-neighbor satisfying cond, call
  // updateAtomic.
  template <class V, class F, class G>
  inline void decodeOutNgh(V* v, long i, F &f, G &g) {
    long d = v->getOutDegree();
    if (d==0) return;
    granular_for(j, 0, d-1, (d-1 > 1000), {
      uintE ngh = v->getOutNeighbor(j);
      if (f.cond(ngh)) {
      auto m = f.updateAtomic(i,ngh);  //由用户提供的原子操作f
        g(ngh, m); // 生成返回结果--下一个边界
      }
    });
    uintE ngh = v->getOutNeighbor(d-1);
    if (v->CliqueNode==NULL) { //判断是否是团结点
      if (f.cond(ngh)) {
        auto m = f.updateAtomic(i,ngh);
        g(ngh, m);
      }
    } else { //遍历团结点包含的所有邻居
      parallel_for(uintE ki = v->CliqueNode->neighbors[0]; ki < v->CliqueNode->neighbors[1]; ki++) {
        if (ki!= i && f.cond(ki)) {
          auto m = f.updateAtomic(i, ki);
          g(ki, m);
        }
      }
    }
  }

 
  // Used by edgeMapSparse. For each out-neighbor satisfying cond, call
  // updateAtomic.
  template <class V, class F, class G>
  inline void decodeOutNghSparse(V* v, long i, uintT o, F &f, G &g) {
    long d = v->getOutDegree();
    if (d==0) return;
    granular_for(j, 0, d-1, (d-1 > 1000), {
      uintE ngh = v->getOutNeighbor(j);
      if (f.cond(ngh)) {
        auto m = f.updateAtomic(i, ngh);
        g(ngh, o+j, m);
      } else {
        g(ngh, o+j);
      }
    });
    uintE ngh = v->getOutNeighbor(d-1);
    if (v->CliqueNode==NULL) {
      if (f.cond(ngh)) {
        auto m = f.updateAtomic(i, ngh);
        g(ngh, o+d-1, m);
      } else {
        g(ngh, o+d-1);
      }
    } else {
      if (ngh<v->CliqueNode->Tdegree) cout<<"error x."<<endl;
      for (long ki = v->CliqueNode->neighbors[0]; ki < v->CliqueNode->neighbors[1]; ki++) {
        uintE index = ki - v->CliqueNode->neighbors[0];
        if (ki == i) continue;
        if (ki > i) index --;
        if (f.cond(ki)) {
          auto m = f.updateAtomic(i, ki);
          g(ki, o+d-1 + index, m);
          index++;
        } else {
          g(ki, o+d-1 + index);
          index++;
        }
      }
    }
  }

  // Used by edgeMapSparse_no_filter. Sequentially decode the out-neighbors,
  // and compactly write all neighbors satisfying g().
  template <class V, class F, class G>
  inline size_t decodeOutNghSparseSeq(V* v, long i, uintT o, F &f, G &g) {
    long d = v->getOutDegree();
    if (d==0) return 0;
    size_t k = 0;
    for (size_t j=0; j<d-1; j++) {
      uintE ngh = v->getOutNeighbor(j);
      if (f.cond(ngh)) {
        auto m = f.updateAtomic(i, ngh);
        bool wrote = g(ngh, o+k, m);
        if (wrote) { k++; }
      }
    }
    uintE ngh = v->getOutNeighbor(d-1);
    if (v->CliqueNode == NULL) {
      if (f.cond(ngh)) {
        auto m = f.updateAtomic(i, ngh);
        bool wrote = g(ngh, o+k, m);
        if (wrote) { k++; }
      }
    } else {
      if (ngh<v->CliqueNode->Tdegree) cout<<"error x."<<endl;
      for (uintE ki = v->CliqueNode->neighbors[0]; ki < v->CliqueNode->neighbors[1]; ki++) {
        if (ki!=i && f.cond(ki)) {
          auto m = f.updateAtomic(i, ki);
          bool wrote = g(ki, o+k, m);
          if (wrote) {k++;}
        }
      }
    }
    return k;
  }


  // Decode the out-neighbors of v, and return the number of neighbors
  // that satisfy f.
  template <class V, class F>
  inline size_t countOutNgh(V* v, long vtx_id, F& f) {
    long d = v->getOutDegree();
    if (d==0) return 0;
    if (d < 2000) {
      size_t ct = 0;
      for (size_t i=0; i<d-1; i++) {
        uintE ngh = v->getOutNeighbor(i);
        if (f(vtx_id, ngh))
          ct++;
      }
      uintE ngh = v->getOutNeighbor(d-1);
      if (v->CliqueNode==NULL) {
        if (f(vtx_id, ngh))
          ct++;
      } else {
        if (ngh<v->CliqueNode->Tdegree) cout<<"error x."<<endl;
        for (uintE ki = v->CliqueNode->neighbors[0]; ki < v->CliqueNode->neighbors[1]; ki++) {
          if (ki!=vtx_id && f(vtx_id, ki)) 
            ct++;
        }
      }
      return ct;
    } else {
      size_t b_size = 2000;
      size_t blocks = 1 + ((d-1)/b_size);
      auto cts = array_imap<uintE>(blocks, [&] (size_t i) { return 0; });
      parallel_for_1(size_t i=0; i<blocks-1; i++) {
        size_t s = b_size*i;
        size_t e = std::min(s + b_size, (size_t)d);
        uintE ct = 0;
        for (size_t j = s; j < e; j++) {
          uintE ngh = v->getOutNeighbor(j);
          if (f(vtx_id, ngh))
            ct++;
        }
        cts[i] = ct;
      }
      size_t s = b_size*(blocks-1);
      size_t e = (size_t)d;
      uintE ct = 0;
      for (size_t j = s; j < e-1; j++) {
        uintE ngh = v->getOutNeighbor(j);
        if (f(vtx_id, ngh))
          ct++;
      }
      uintE ngh = v->getOutNeighbor(e-1);
      if (v->CliqueNode == NULL) {
        if (f(vtx_id, ngh))
          ct++;
      } else {
        if (ngh<v->CliqueNode->Tdegree) cout<<"error x."<<endl;
        for (uintE ki = v->CliqueNode->neighbors[0]; ki < v->CliqueNode->neighbors[1]; ki++) {
          if (f(vtx_id, ki))
            ct++;
        }
      }
      cts[blocks-1] = ct;
      size_t count = 0;
      return pbbs::reduce_add(cts);
    }
  }



  // Decode the out-neighbors of v. Apply f(src, ngh) and store the result
  // using g.
  template <class V, class E, class F, class G>
  inline void copyOutNgh(V* v, long src, uintT o, F& f, G& g) {
    long d = v->getOutDegree();
    if (d==0) return;
    granular_for(j, 0, d-1, (d > 1000), {
      uintE ngh = v->getOutNeighbor(j);
      E val = f(src, ngh);
      g(ngh, o+j, val);
    });
    uintE ngh = v->getOutNeighbor(d-1);
    if (v->CliqueNode==NULL) {
      E val = f(src, ngh);
      g(ngh, o+d-1, val);
    } else {
      if (ngh<v->CliqueNode->Tdegree) cout<<"error x."<<endl;
      long index = 0;
      for(uintE ki = v->CliqueNode->neighbors[0]; ki < v->CliqueNode->neighbors[1]; ki++) {
        if (ki == src) continue;
        E val = f(src, ki);
        g(ki, o+d-1+index, val);
        index++;
      }
    }
  }


  // TODO(laxmand): Add support for weighted graphs.
  template <class V, class Pred>
  inline size_t packOutNgh(V* v, long vtx_id, Pred& p, bool* bits, uintE* tmp) {
    cout<<"-forbid-"<<endl;
    uintE d = v->getOutDegree();
    if (d < 5000) {
      size_t k = 0;
      for (size_t i=0; i<d; i++) {
        uintE ngh = v->getOutNeighbor(i);
        if (p(vtx_id, ngh)) {
          v->setOutNeighbor(k, ngh);
          k++;
        }
      }
      v->setOutDegree(k);
      return k;
    } else {
      parallel_for(size_t i=0; i<d; i++) {
        uintE ngh = v->getOutNeighbor(i);
        tmp[i] = ngh;
        bits[i] = p(vtx_id, ngh);
      }
      size_t k = sequence::pack(tmp, v->getOutNeighbors(), bits, d);
      v->setOutDegree(k);
      return k;
    }
  }
}

struct zipSymmetricVertex {
  uintE* neighbors;
  uintT degree;
  zipSymmetricVertex *CliqueNode;
  uintT Tdegree;
  void del() {
    cout<<"hello.."<<endl;
    // if (neighbors != NULL)
    // free(neighbors); 
 }

zipSymmetricVertex(uintE* n, uintT d, uintT t)
: neighbors(n), degree(d) , Tdegree(t) {}

  uintE* getInNeighbors () { return neighbors; }
  const uintE* getInNeighbors () const { return neighbors; }
  uintE* getOutNeighbors () { return neighbors; }
  const uintE* getOutNeighbors () const { return neighbors; }
  uintE getInNeighbor(uintT j) const { if (j >= degree) cout<<"error!!!!"<<endl; return neighbors[j]; }
  uintE getOutNeighbor(uintT j) const { if (j >= degree) cout<<"error!!!!"<<endl; return neighbors[j]; }

  void setInNeighbor(uintT j, uintE ngh) { neighbors[j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { neighbors[j] = ngh; }
  void setInNeighbors(uintE* _i) { neighbors = _i; }
  void setOutNeighbors(uintE* _i) { neighbors = _i; }


  uintT getInDegree() const { return degree; }
  uintT getOutDegree() const { return degree; }
  uintT getTrueInDegree() const { return Tdegree; }
  uintT getTrueOutDegree() const { return Tdegree; }
  void setInDegree(uintT _d) { degree = Tdegree = _d; }
  void setOutDegree(uintT _d) { degree = Tdegree = _d; }
  void flipEdges() {}

  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(long v_id, VS& vertexSubset, F &f, G &g, bool parallel = 0) {
    decode_zip::decodeInNghBreakEarly<zipSymmetricVertex, F, G, VS>(this, v_id, vertexSubset, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(long i, F &f, G& g) {
     decode_zip::decodeOutNgh<zipSymmetricVertex, F, G>(this, i, f, g);
  }

  template <class F, class G>
  inline void decodeOutNghSparse(long i, uintT o, F &f, G &g) {
    decode_zip::decodeOutNghSparse<zipSymmetricVertex, F>(this, i, o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(long i, uintT o, F &f, G &g) {
    return decode_zip::decodeOutNghSparseSeq<zipSymmetricVertex, F>(this, i, o, f, g);
  }

  template <class E, class F, class G>
  inline void copyOutNgh(long i, uintT o, F& f, G& g) {
    decode_zip::copyOutNgh<zipSymmetricVertex, E>(this, i, o, f, g);
  }

  template <class F>
  inline size_t countOutNgh(long i, F &f) {
    return decode_zip::countOutNgh<zipSymmetricVertex, F>(this, i, f);
  }

  template <class F>
  inline size_t packOutNgh(long i, F &f, bool* bits, uintE* tmp1, uintE* tmp2) {
    return decode_zip::packOutNgh<zipSymmetricVertex, F>(this, i, f, bits, tmp1);
  }

};

struct zipAsymmetricVertex {
  uintE* neighbors;
#ifndef WEIGHTED
  uintE* inNeighbors, *outNeighbors;
#else
  intE* inNeighbors, *outNeighbors;
#endif
  uintT outDegree;
  uintT inDegree;
  uintT degree;
  zipSymmetricVertex *CliqueNode;
  uintT Tdegree;
  void del() {free(inNeighbors); free(outNeighbors);}
#ifndef WEIGHTED
zipAsymmetricVertex(uintE* iN, uintE* oN, uintT id, uintT od)
#else
zipAsymmetricVertex(intE* iN, intE* oN, uintT id, uintT od)
#endif
: inNeighbors(iN), outNeighbors(oN), inDegree(id), outDegree(od) {}
#ifndef WEIGHTED
  uintE* getInNeighbors () { return inNeighbors; }
  const uintE* getInNeighbors () const { return inNeighbors; }
  uintE* getOutNeighbors () { return outNeighbors; }
  const uintE* getOutNeighbors () const { return outNeighbors; }
  uintE getInNeighbor(uintT j) const { return inNeighbors[j]; }
  uintE getOutNeighbor(uintT j) const { return outNeighbors[j]; }
  void setInNeighbor(uintT j, uintE ngh) { inNeighbors[j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[j] = ngh; }
  void setInNeighbors(uintE* _i) { inNeighbors = _i; }
  void setOutNeighbors(uintE* _i) { outNeighbors = _i; }
#else
  intE* getInNeighbors () { return inNeighbors; }
  const intE* getInNeighbors () const { return inNeighbors; }
  intE* getOutNeighbors () { return outNeighbors; }
  const intE* getOutNeighbors () const { return outNeighbors; }
  intE getInNeighbor(uintT j) const { return inNeighbors[2*j]; }
  intE getOutNeighbor(uintT j) const { return outNeighbors[2*j]; }
  intE getInWeight(uintT j) const { return inNeighbors[2*j+1]; }
  intE getOutWeight(uintT j) const { return outNeighbors[2*j+1]; }
  void setInNeighbor(uintT j, uintE ngh) { inNeighbors[2*j] = ngh; }
  void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[2*j] = ngh; }
  void setInWeight(uintT j, uintE wgh) { inNeighbors[2*j+1] = wgh; }
  void setOutWeight(uintT j, uintE wgh) { outNeighbors[2*j+1] = wgh; }
  void setInNeighbors(intE* _i) { inNeighbors = _i; }
  void setOutNeighbors(intE* _i) { outNeighbors = _i; }
#endif

  uintT getInDegree() const { return inDegree; }
  uintT getOutDegree() const { return outDegree; }
  uintT getTrueInDegree() const { return inDegree; }
  uintT getTrueOutDegree() const { return outDegree; }
  void setInDegree(uintT _d) { inDegree = _d; }
  void setOutDegree(uintT _d) { outDegree = _d; }
  void flipEdges() { swap(inNeighbors,outNeighbors); swap(inDegree,outDegree); }

  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(long v_id, VS& vertexSubset, F &f, G &g, bool parallel = 0) {
    decode_zip::decodeInNghBreakEarly<zipAsymmetricVertex, F, G, VS>(this, v_id, vertexSubset, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(long i, F &f, G &g) {
    decode_zip::decodeOutNgh<zipAsymmetricVertex, F, G>(this, i, f, g);
  }

  template <class F, class G>
  inline void decodeOutNghSparse(long i, uintT o, F &f, G &g) {
    decode_zip::decodeOutNghSparse<zipAsymmetricVertex, F>(this, i, o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(long i, uintT o, F &f, G &g) {
    return decode_zip::decodeOutNghSparseSeq<zipAsymmetricVertex, F>(this, i, o, f, g);
  }

  template <class E, class F, class G>
  inline void copyOutNgh(long i, uintT o, F& f, G& g) {
    decode_zip::copyOutNgh<zipAsymmetricVertex, E>(this, i, o, f, g);
  }

  template <class F>
  inline size_t countOutNgh(long i, F &f) {
    return decode_zip::countOutNgh<zipAsymmetricVertex, F>(this, i, f);
  }

  template <class F>
  inline size_t packOutNgh(long i, F &f, bool* bits, uintE* tmp1, uintE* tmp2) {
    return decode_zip::packOutNgh<zipAsymmetricVertex, F>(this, i, f, bits, tmp1);
  }

};

#endif
