
#ifndef MESHCLASS_H
#define MESHCLASS_H

#include <vector>
#include <vcg/complex/complex.h>

class MyFace;
class MyVertex;

struct MyUsedTypes : public vcg::UsedTypes< vcg::Use<MyVertex>::AsVertexType, vcg::Use<MyFace>::AsFaceType> {
};

class MyVertex : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::BitFlags, vcg::vertex::Qualityf > { //vcg::vertex::Normal3f, 
};

class MyFace : public vcg::Face < MyUsedTypes, vcg::face::VertexRef, vcg::face::FFAdj, vcg::face::BitFlags,  vcg::face::Qualityf > { //vcg::face::Normal3f, 
};

class MyMesh : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {
};

typedef MyMesh MeshType;
typedef typename MeshType::ScalarType ScalarType;
typedef typename MeshType::CoordType CoordType;
typedef typename MeshType::FaceType FaceType;
typedef typename MeshType::FaceType* FaceTypeP;
typedef typename MeshType::VertexType VertexType;
typedef typename MeshType::FaceIterator FaceIterator;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::ConstFaceIterator ConstFaceIterator;

#endif /* MESHCLASS_H */

