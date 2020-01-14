

#include <iostream>
#include <fstream>
#include "help_functions_for_tests.h"
//#include "ml_mesh_type.h"
using std::cout;
using std::endl;
using std::string;

int main(int argc, char**argv) {
  string meshname(argv[1]);
  cout << "Abrindo " << meshname << "..." << endl;
    
  MeshType m, m2;
  int err = vcg::tri::io::ImporterOFF<MeshType>::Open(m, meshname.c_str()); //, mask );
  if (err) {
      std::cerr << "Unable to open mesh " << meshname << " : " << vcg::tri::io::ImporterOFF<MeshType>::ErrorMsg(err) << std::endl;
      exit(-1);
  }
  
  cout << "Malha carregada. " << endl;
  cout << "  Nº faces=" << m.FN() << endl;
  cout << "  Nº vertices=" << m.VN() << endl;
  cout << "  Nº MemUsed=" << m.MemUsed() << " bytes" << endl; //unidade do sizeof
 
  vcg::tri::UpdateTopology<MeshType>::FaceFace(m); //atualizar ponteiros para faces vizinhas e numeração das arestas
  
  //checks for consistency
  if (!vcg::tri::Clean<MeshType>::IsBitTriOnly(m))
  {
    cout << "Malha com arestas falsas. Selecione uma malha tri-only"<<endl;
    exit(-2);
  }
  if (vcg::tri::Clean<MeshType>::CountNonManifoldVertexFF(m,false) != 0) 
  {
    cout << "Malha com vertex non manifold"<<endl;
    exit(-3);
  }
  if (vcg::tri::Clean<MeshType>::CountNonManifoldEdgeFF(m,false) != 0)
  {
    cout << "Malha com edge non manifold"<<endl;
    exit(-4);
  }
  
  vcg::tri::Append<MeshType, MeshType>::MeshCopy(m2, m, false, true); //salvar copia em m2
  std::chrono::high_resolution_clock::time_point start; //tempo antes da função. It is the clock type with the highest precision.
  std::chrono::high_resolution_clock::time_point stop; //tempo depois da função
  int nclusters=0; //número de clusters retornado pelas funções
  int nIter=10; //number of iterations
  vector<vcg::Histogram<ScalarType>> histFaces,histVerts; //histograms of faces quality and vertex valence
  int qualityComp = 0; //qualidade das faces para comparação. distância para 90º. A qualidade dos vértices é a valência.
  vector<long int> tempos; //duração de cada função de agrupamento
  vector<int> clusters; //
  vector<int> facesFinal;
  vector<int> vertsFinal;
  //int metodo;
  
  int nAlg=10; //nº de algoritmos
  string nomes[nAlg]={
                  "Algoritmo 0", //"No clusters" 
                  "Algoritmo 1", //"Maior Qualidade Angular", 
                  "Algoritmo 2", //"Maior aresta (Luis Velho)", 
                  "Algoritmo 3", //"Maior Qualidade Orto Planar", 
                  "Algoritmo 4", //"Menor Tri Aspect Ratio Maior Aresta", 
                  "Algoritmo 5", //"Menor Tri Aspect Ratio Maior Aresta com TestConvex", 
                  "Algoritmo 6", //"Linear Level 0", 
                  "Algoritmo 7", //"Linear Level 1", 
                  "Algoritmo 8", //"Linear Level 2", 
                  "Algoritmo 9", //"Linear Level 3", 
                  };
  
  assert(vcg::tri::BitQuad<MeshType>::CountFauxs(m)==0);
  
  for (int metodo=0; metodo<nAlg; metodo++)
  {
    cout <<"**Executando "<<nomes[metodo]<<"..."<<endl;
    
    for(int i=0;i<nIter;++i) //iterações
    {
      vcg::tri::Append<MeshType, MeshType>::MeshCopy(m, m2, false, true); //restaura malha inicial
      assert(vcg::tri::BitQuad<MeshType>::CountFauxs(m)==0); //apagar
      start = std::chrono::high_resolution_clock::now();
        nclusters = ApplyPairing(m, metodo); //aplica metodos de agrupamento
      stop = std::chrono::high_resolution_clock::now();
      tempos.push_back(std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()); //salva duração
      
      if (metodo == 6 || metodo == 7 || metodo == 8 || metodo == 9){
        clusters.push_back(vcg::tri::BitQuad<MeshType>::CountFauxs(m));
        assert(clusters.back() != 0);
      }else{
        assert(nclusters == vcg::tri::BitQuad<MeshType>::CountFauxs(m));
        clusters.push_back(nclusters);
      }
    }
    
    SaveMesh(m, nomes[metodo] + " - agrupamentos"); //salva malha após criação de clusters
    GenerateHistogramsData(m,nomes[metodo], qualityComp, histFaces, histVerts); //aplica refinamento e salva qualidades em arquivos
    facesFinal.push_back(m.FN()); //salva número de faces na malha final
    vertsFinal.push_back(m.VN()); //salva número de vértices na malha final
    SaveMesh(m, nomes[metodo] + " - final"); //salvar malha após refinamento

  }
  
  cout<<"\nQtd de iterações="<<tempos.size()<<endl;
  

  for (int i = 0; i < nAlg; i++) { //ordena as medidas de cada algorítimo
    std::sort(tempos.begin()+i*nIter, tempos.begin()+(i+1)*nIter);
  }
//
//  for (auto x : tempos) {
//    cout<<x<<" ";
//  }
//  cout<<endl<<"clusters: ";
//  for (auto x : clusters) {
//    cout<<x<<" ";
//  }
//  
  cout<<endl<<"Resultados: "<<endl;
  if (nIter%2==0){
    int p=0;
    for (int i = 0; i < nAlg; i++) {
      p=nIter/2+i*nIter;
      cout<<nomes[i]<< endl;
      //cout<<"  Faces triangulares iniciais="<< m2.FN()<<endl;
      //cout<<"  Vértices iniciais="<< m2.VN()<<endl;
      cout<<"  Agrupamentos realizados="<<clusters[p]<<endl;
      cout<<"  Vértices finais="<<vertsFinal[i]<<endl; //
      cout<<"  Quadriláteros finais="<<facesFinal[i]/2<<endl; //assert(facesFinal[i]/2 == res[p]*4+(m2.FN()-res[p]*2)*3);//
      cout<<"  Mediana dos tempos de execução do agrupamento(micros)="<< (tempos[p]+tempos[p-1])/2<<endl;
      cout<<"  Média das qualidades das faces="<<histFaces[i].Avg()<<endl;
      cout<<"  Desvio padrão das qualidades das faces=" << histFaces[i].StandardDeviation()<<endl;
      cout<<"  Qualidade da Pior face="<<histFaces[i].MinElem() << endl;
      //cout<<"  Percentil(0.5)="<<histFaces[i].Percentile(0.5)<<endl;
      cout<<"  Quantidade de vértices regulares=" << histVerts[i].BinCount(4)<<"("<<histVerts[i].BinCount(4)/histVerts[i].Cnt()<<")"<<endl;
      //if (histVerts[i].BinCountInd(9) > 0) cout<< "*bin9*="<<histVerts[i].BinCountInd(9)<<endl;
      //if (histVerts[i].BinCountInd(10) > 0) cout<< "*bin10*="<<histVerts[i].BinCountInd(10)<<endl;
    }
  }else{
//    int p=0;
//    for (int i = 0; i < nAlg; i++) {
//      p=nIter/2+i*nIter;
//      cout<<nomes[i]<< ": clusters="<<clusters[p]<<"  Mediana (micros)="<< tempos[p]<<endl;
//    }
  }
  cout<<endl;
  
  //cout<<"nº histograma faces: "<<histFaces.size()<<endl;
  //cout<<"nº histograma vertex: "<<histVerts.size()<<endl;
  //histVerts[0].FileWrite("zv");
  //histFaces[0].FileWrite("zf");
  
//  for (int i=0; i<8; i++){
//    cout<<histFaces[i].Percentile(0.5)<< " ";
//  }
  //media
  //desvio padrão
  //.Percentile(0.5)
  
//  vcg::tri::Append<MeshType, MeshType>::MeshCopy(m, m2, false, true);
//  vcg::tri::UpdateQuality<MeshType>::VertexValence(m);
//  vcg::Histogram<ScalarType> Hv;
//  ComputePerVertexQualityHistogram2(m,Hv);
//  Hv.FileWrite("histograma vertex tri only");
//  if (Hv.BinCountInd(10) > 0) cout<< "\n\n*bin10 tri only*="<<Hv.BinCountInd(10)<<endl;

  return 0;
}
