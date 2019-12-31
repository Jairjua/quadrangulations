#include "help_functions_for_tests.h"
#include "eigenlib/Eigen/src/StlSupport/StdVector.h"
//#include <string>
//#include <tuple>
//#include <vector>
//#include <cmath>
//#include <utility>      // std::pair
//#include <chrono>       //high_resolution_clock::now()
//#include <fstream>
//
//
//#include <wrap/io_trimesh/import.h>

//#include <wrap/io_trimesh/export_obj.h>

//
//#define PI 3.1415926535897
//
using std::cout;
using std::endl;
using std::string;
//typedef vcg::tri::GeometricInterpolator<typename MeshType::VertexType> Interpolator;
//typedef vcg::tri::BitQuad<MeshType> BQ; 

int ApplyPairing(MeshType& m, int nalg)
{
  switch (nalg)
  {
    case 0:
      return 0;
      break;
    case 1:
      return vcg::tri::BitQuadCreation<MeshType>::MakeDominantByOrder(m,0); //maior qualidade orto planar
      break;
    case 2:
      return vcg::tri::BitQuadCreation<MeshType>::MakeDominantByOrder(m,1); //Luis Velho
      break;
    case 3:
      return vcg::tri::BitQuadCreation<MeshType>::MakeDominantByOrder(m,2); //Maior Qualidade Orto Planar
      break;  
    case 4:
      return vcg::tri::BitQuadCreation<MeshType>::MakeDominantByTriOrder(m,1); //Menor Tri Aspect Ratio Maior Aresta
      break;
    case 5:
      return vcg::tri::BitQuadCreation<MeshType>::MakeDominantByTriOrder(m,1,0,false,true); //Menor Tri Aspect Ratio Maior Aresta com TestConvex
      break;
    case 6:
      vcg::tri::BitQuadCreation<MeshType>::MakeDominant(m,0); //linear level 0 - 1 iteração
      return -1;
      break;
    case 7:
      vcg::tri::BitQuadCreation<MeshType>::MakeDominant(m,2); //linear level 2 - 4 iterações
      return -1;
      break;
    case 8:
      vcg::tri::BitQuadCreation<MeshType>::MakeDominant(m,1); //linear level 1 - 3 iterações
      return -1;
      break;
    case 9:
      vcg::tri::BitQuadCreation<MeshType>::MakeDominant(m,3); //linear level 3 - 2 iterações
      return -1;
      break;
    default:
      break;
  }
}

void GenerateHistogramsData (MeshType& m, string identifier, int faceQuality, 
                             vector<vcg::Histogram<ScalarType>>& histFaces, vector<vcg::Histogram<ScalarType>>& histVerts)
{
  //cout<<"faces antes refinamento=" <<m.FN() <<endl;
  vcg::tri::BitQuadCreation<MeshType>::MakePureByCatmullClark(m); //refinamento 
  //cout <<"faces após refinamento="<<m.FN()<< " vertex depois="<<m.VN()<<endl;
  
  //cout <<"valency quality antes="<<m.vert.begin()->Q()<<endl;
  vcg::tri::BitQuad<MeshType>::UpdateValencyInQualityWithBorders(m);
  //cout <<"valency quality depois="<<m.vert.begin()->Q()<<endl;
  //cout <<"quality antes="<<m.face.begin()->Q()<<endl;
  vcg::tri::BitQuad<MeshType>::MeasureQuality(m,faceQuality); 
  //cout <<"quality depois="<<m.face.begin()->Q()<<endl;
  
  vcg::Histogram<ScalarType> Hf,Hv;
  ComputePerFaceQualityHistogram2(m,Hf,false,20);
  ComputePerVertexQualityHistogram2(m,Hv);
  
  histFaces.push_back(Hf);
  histVerts.push_back(Hv);
  Hf.FileWrite(identifier + " - Histograma Faces");
  Hv.FileWrite(identifier + " - Histograma Vertices");
  //SaveQualitys(m, identifier+" - face qualitys", identifier+" - vertex qualitys");
  
}

 void ComputePerFaceQualityHistogram2( MeshType & m, vcg::Histogram<ScalarType> &h, bool selectionOnly,int HistSize)
{
  vcg::tri::RequirePerFaceQuality(m);
  std::pair<ScalarType, ScalarType> minmax = vcg::tri::Stat<MeshType>::ComputePerFaceQualityMinMax(m);
  h.Clear();
  if (minmax.first <= 0 || minmax.second > 1){
    h.SetRange( minmax.first-0.01, minmax.second+0.1, HistSize );
  }else{
    h.SetRange( 0, 1, HistSize );
  }
  for(FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    if(!(*fi).IsD() &&  ((!selectionOnly) || (*fi).IsS()) ){
      assert(!vcg::math::IsNAN((*fi).Q()) && "You should never try to compute Histogram with Invalid Floating points numbers (NaN)");
      h.Add((*fi).Q());
    }
}

void ComputePerVertexQualityHistogram2( MeshType & m, vcg::Histogram<ScalarType> &h, bool selectionOnly, int HistSize)    // V1.0
{
  vcg::tri::RequirePerVertexQuality(m);
  //std::pair<ScalarType, ScalarType> minmax = ComputePerVertexQualityMinMax(m);

  h.Clear();
  h.SetRange( 0,9, 9);
  for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
    if(!(*vi).IsD() &&  ((!selectionOnly) || (*vi).IsS()) )
    {
      assert(!vcg::math::IsNAN((*vi).Q()) && "You should never try to compute Histogram with Invalid Floating points numbers (NaN)");
      h.Add((*vi).Q());
    }
  // Sanity check; If some very wrong value has happened in the Q value,
  // the histogram is messed. If a significant percentage (20% )of the values are all in a single bin
  // we should try to solve the problem. No easy solution here.
  // We choose to compute the get the 1percentile and 99 percentile values as new mixmax ranges
  // and just to be sure enlarge the Histogram.

//  if(h.MaxCount() > HistSize/5)
//  {
//    std::vector<ScalarType> QV;
//    QV.reserve(m.vn);
//    for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
//      if(!(*vi).IsD()) QV.push_back((*vi).Q());
//
//    std::nth_element(QV.begin(),QV.begin()+m.vn/100,QV.end());
//    ScalarType newmin=*(QV.begin()+m.vn/100);
//    std::nth_element(QV.begin(),QV.begin()+m.vn-m.vn/100,QV.end());
//    ScalarType newmax=*(QV.begin()+m.vn-m.vn/100);
//
//    h.Clear();
//    h.SetRange(newmin, newmax, HistSize*50);
//    for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
//      if(!(*vi).IsD() && ((!selectionOnly) || (*vi).IsS()) )
//        h.Add((*vi).Q());
//  }
}
 
void SaveQualitys(MeshType& m, string nameFileFaces, string nameFileVert)
{

  std::ofstream fileFace, fileVert;
  std::vector<ScalarType> vf, vv;

  for (FaceIterator fi = m.face.begin();  fi!=m.face.end(); fi++){
    vf.push_back(fi->Q());
  }
  sort(vf.begin(),vf.end());
  
  for (VertexIterator vi = m.vert.begin();  vi!=m.vert.end(); vi++){
    vv.push_back(vi->Q());
  }
  sort(vv.begin(),vv.end());
  
  fileFace.open(nameFileFaces);//, ios::app);
  if(!fileFace.is_open()){ //testar deixar aberto sempre
      cout<<"não deu pra abrir o arquivo: "<< nameFileFaces<<endl;
      return;
  }
  //save face qualitys
  for (int i = 0; i < vf.size(); i++)
  {
    //cout << vf[i]<< " ";
    if (i%2==0){
      fileFace << vf[i] << "\n";
      assert (vf[i]-vf[i+1] < 0.001 && vf[i+1] - vf[i] < 0.001); //entradas são suplicadas. a igualdade n funciona provavelmente pq o q é calculado para cada fac e não copiado
    }
  }
  fileFace.close();
  
  cout<< endl;
  fileVert.open(nameFileVert);//, ios::app);
  if(!fileVert.is_open()){ //testar deixar aberto sempre
      cout<<"não deu pra abrir o arquivo: "<<nameFileVert<<endl;
      return;
  }
  //save vertex qualitys
  for (int i = 0; i < vv.size(); i++)
  {
    //cout << vv[i]<< " ";
    fileVert << vv[i] << "\n";
  }
  fileVert.close();
}

bool Pair(FaceType* f, int edge)
{
  if (f->FFp(edge) == f) return false;
  f->SetF(edge);
  f->FFp(edge)->SetF(f->FFi(edge));
  return true;
}

int CountTriAndQuad(MeshType& m, int& nQuads, int& nTris)
{
  if(!vcg::tri::BitQuadCreation<MeshType>::IsTriQuadOnly(m)){
    return -2;
  }
  nQuads = vcg::tri::BitQuad<MeshType>::CountFauxs(m);
  if (nQuads >=0) nTris = m.FN() - 2 * nQuads;
  return nQuads;
}

void PrintVertexCoord(const MeshType& m)
{//MeshType::VertexIterator
  int n=0;
  for (auto vi = m.vert.begin(); vi != m.vert.end(); vi++)
  {
    cout << "vert-" << n << " \t";
    for (int i = 0; i < 3; ++i)
    {
      cout << vi->P()[i] << ","; //colocar numeração com 0 antes.
    }
    n++;
    cout << "\b " << endl; // apagar a última vírgula antes de pular
  }

}

string CreateMesh(MeshType& m, int i)
{
  string text[10] ={"Lista de descrição das malhas",\
    "Mesh tipo 1(Lista de vértices não convexo:4v2f)", \
    "Mesh tipo 2(Lista de vértices não convexo e impar:5v3f)",\
    "Mesh tipo 3(Lista de vértices:7v6f)",\
    "Mesh tipo 4(Lista de vértices sem borda:6v8f)",\
    "Mesh tipo 5(Tetrahedron)",\
    "Mesh tipo 6(Octahedron)"    
  };

  switch (i) {
  case 0:
  {
    return string("No mesh created");
    
    break;
  }
  case 1:
  {
    std::vector<vcg::Point3f> coordVec;
    std::vector<vcg::Point3i> indexVec;
    coordVec.push_back(vcg::Point3f(0, 0, 0));
    coordVec.push_back(vcg::Point3f(2, 0, 0));
    coordVec.push_back(vcg::Point3f(2, 2, 0));
    coordVec.push_back(vcg::Point3f(1, 1, 0));

    indexVec.push_back(vcg::Point3i(0, 1, 3));
    indexVec.push_back(vcg::Point3i(1, 2, 3));

    vcg::tri::BuildMeshFromCoordVectorIndexVector(m, coordVec, indexVec);

    break;
  }
  case 2:
  {
    std::vector<vcg::Point3f> coordVec;
    std::vector<vcg::Point3i> indexVec;
    coordVec.push_back(vcg::Point3f(0, 0, 0));
    coordVec.push_back(vcg::Point3f(1, 0, 0));
    coordVec.push_back(vcg::Point3f(1.5, 2, 0));
    coordVec.push_back(vcg::Point3f(1, 1, 0));
    coordVec.push_back(vcg::Point3f(0, 2, 0));

    indexVec.push_back(vcg::Point3i(0, 1, 3));
    indexVec.push_back(vcg::Point3i(1, 2, 3));
    indexVec.push_back(vcg::Point3i(0, 3, 4));

    vcg::tri::BuildMeshFromCoordVectorIndexVector(m, coordVec, indexVec);
    
    break;
  }
  case 3:
  {
    std::vector<vcg::Point3f> coordVec;
    std::vector<vcg::Point3i> indexVec;
    coordVec.push_back(vcg::Point3f(0, 0, 0));
    coordVec.push_back(vcg::Point3f(1, 0, 0));
    coordVec.push_back(vcg::Point3f(1.5, 0.5, 0));
    coordVec.push_back(vcg::Point3f(1, 1, 0));
    coordVec.push_back(vcg::Point3f(0, 1, 0));
    coordVec.push_back(vcg::Point3f(0.25, 1.5, 0));
    coordVec.push_back(vcg::Point3f(2, 1.75, 0));

    indexVec.push_back(vcg::Point3i(0, 1, 3));
    indexVec.push_back(vcg::Point3i(1, 2, 3));
    indexVec.push_back(vcg::Point3i(0, 3, 4));
    indexVec.push_back(vcg::Point3i(4, 3, 5));
    indexVec.push_back(vcg::Point3i(3, 6, 5));
    indexVec.push_back(vcg::Point3i(3, 2, 6));

    vcg::tri::BuildMeshFromCoordVectorIndexVector(m, coordVec, indexVec);

    break;
  }
  case 4:
  {
    std::vector<vcg::Point3f> coordVec;
    std::vector<vcg::Point3i> indexVec;
    coordVec.push_back(vcg::Point3f(0, 0, 0));
    coordVec.push_back(vcg::Point3f(1, 0, 0));
    coordVec.push_back(vcg::Point3f(1, 1, 0));
    coordVec.push_back(vcg::Point3f(1, 0, -1));
    coordVec.push_back(vcg::Point3f(0, 1, -1));
    coordVec.push_back(vcg::Point3f(0, 0, -1));

    indexVec.push_back(vcg::Point3i(0, 1, 2));
    indexVec.push_back(vcg::Point3i(1, 3, 2));
    indexVec.push_back(vcg::Point3i(3, 5, 2));
    indexVec.push_back(vcg::Point3i(5, 4, 2));
    indexVec.push_back(vcg::Point3i(4, 0, 2));

    indexVec.push_back(vcg::Point3i(0, 4, 5));
    indexVec.push_back(vcg::Point3i(0, 5, 1));
    indexVec.push_back(vcg::Point3i(1, 5, 3));

    vcg::tri::BuildMeshFromCoordVectorIndexVector(m, coordVec, indexVec);

    break;
  }
  case 5:
  {
    vcg::tri::Tetrahedron(m);
    
    break;
  }
  case 6:
  {
    vcg::tri::Octahedron(m);
    
    break;
  }
  default:
    break;
  }
  return text[i];
}

string CreateTestMesh(MeshType& m,
                     int nMesh,
                     bool showCoord, //and descrição
                     bool updateFFAdj,
                     bool updateNormalFace)
{
  string desc = CreateMesh(m,nMesh);
  
  if (showCoord) 
  {
    cout << desc << endl;
    PrintVertexCoord(m);
  }
  if (updateFFAdj)  vcg::tri::UpdateTopology<MeshType>::FaceFace(m);
  //if (updateNormalFace) vcg::tri::UpdateNormal<MeshType>::PerFace(m);
  
  return desc;
}

void SaveMesh(MeshType& m, string name, int type){ //type=4
    switch (type) {
        case 1: //save mesh in off format with Mask::IOM_BITPOLYGONAL
            vcg::tri::io::ExporterOFF<MeshType>::Save(m, name.append(".off").c_str(),vcg::tri::io::Mask::IOM_BITPOLYGONAL);
            break;
        case 2: //save trimesh in off format
            vcg::tri::io::ExporterOFF<MeshType>::Save(m, name.append(".off").c_str());
            break;
        case 3: //save in VMI dump format
            vcg::tri::io::ExporterVMI<MeshType>::Save(m, name.append(".vmi").c_str());
//            int mask = vcg::tri::io::Mask::IOM_FACEFLAGS | vcg::tri::io::Mask::IOM_FACEINDEX | vcg::tri::io::Mask::IOM_VERTCOORD;
//            vcg::tri::io::ImporterVMI<MeshType>::Open(m2, "file.vmi", mask);
            break;
        case 4: //save mesh in off format with vertex in order with Mask::IOM_BITPOLYGONAL
            vcg::tri::io::ExporterOFF<MeshType>::Save(m, name.append(".off").c_str(),vcg::tri::io::Mask::IOM_BITPOLYGONAL, true);
            break;
        default:
            break;
    }
}
