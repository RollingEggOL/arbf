//
// Created by Ke Liu on 8/3/16.
//

#ifndef ARBF_TEST_MESHFACTORY_H
#define ARBF_TEST_MESHFACTORY_H

#include <memory>
#include "Mesh.h"

/*
 * Abstract factory to create mesh.
 */
class MeshFactory {
public:
    virtual ~MeshFactory();
    virtual std::unique_ptr<Mesh> createMeshFromFile(const char* filename) = 0;
};

/*
 * Concrete factory to create mesh.
 */
class TriMeshFactory: public MeshFactory {
public:
    ~TriMeshFactory();
    virtual std::unique_ptr<Mesh> createMeshFromFile(const char* filename);
};

/*
 * Abstract mesh factory decorator.
 */
class MeshFactoryWithExperimentalFeatures: public MeshFactory {
public:
    virtual ~MeshFactoryWithExperimentalFeatures();
    virtual std::unique_ptr<Mesh> createMeshFromFile(const char* filename);
    void setMeshFactory(std::unique_ptr<MeshFactory> factory);

private:
    std::unique_ptr<MeshFactory> m_factory;
};

/*
 * Concrete mesh factory decorator
 */
class MeshFactoryWithEnlargedBoundingBox: public MeshFactoryWithExperimentalFeatures {
public:
    virtual std::unique_ptr<Mesh> createMeshFromFile(const char* filename);

private:
    void enlargeBoundingBox(Mesh* mesh);
};

#endif //ARBF_TEST_MESHFACTORY_H
