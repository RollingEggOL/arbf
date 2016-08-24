//
// Created by Ke Liu on 8/4/16.
//

#include <vector>
#include <set>
#include <cstdlib>
#include <cstdio>
#include "../include/MeshNeighborhood.h"

// MeshNeighborhood implementations
MeshNeighborhood::~MeshNeighborhood() {}


// OneRingMeshNeighborhood implementations
OneRingMeshNeighborhood::~OneRingMeshNeighborhood() {}

void OneRingMeshNeighborhood::computeVertexNeighbors(const Mesh *mesh) {
    m_vv_1ring.resize(mesh->getNumVertices());
    m_vf_1ring.resize(mesh->getNumVertices());

    // build 1-ring VERTEX <-> FACE mapping
    for (int f = 0; f < mesh->getNumFaces(); ++f) {
        m_vf_1ring[mesh->getFaces()[f].a].push_back(f);
        m_vf_1ring[mesh->getFaces()[f].b].push_back(f);
        m_vf_1ring[mesh->getFaces()[f].c].push_back(f);
    }

    // find 1-ring VERTEX <-> VERTEX mapping
    for (int v = 0; v < mesh->getNumVertices(); ++v) {
        for (size_t j = 0; j < m_vf_1ring[v].size(); ++j) {
            m_vv_1ring[v].insert(mesh->getFaces()[m_vf_1ring[v][j]].a);
            m_vv_1ring[v].insert(mesh->getFaces()[m_vf_1ring[v][j]].b);
            m_vv_1ring[v].insert(mesh->getFaces()[m_vf_1ring[v][j]].c);
        }
    }

//    if (Config::isPrintingDebugInfo) {
//        string filename = string("DEBUG.vf_1ring.txt");
//        FILE *fp = nullptr;
//
//        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
//            fprintf(stderr, "WARNING: open DEBUG.vf_1ring.txt failed ...\n");
//            return;
//        }
//
//        for (int i = 0; i < s_vf_1ring.size(); ++i) {
//            fprintf(fp, "[v <-> f] [%ld] v%d: ", s_vf_1ring[i].size(), i);
//            for (int j = 0; j < s_vf_1ring[i].size(); ++j) {
//                fprintf(fp, "%d ", s_vf_1ring[i][j]);
//            }
//            fprintf(fp, "\n");
//        }
//
//        fclose(fp);
//        fp = nullptr;
//        filename = string("DEBUG.vv_1ring.txt");
//
//        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
//            fprintf(stderr, "WARNING: open DEBUG.vv_1ring.txt failed ...\n");
//            return;
//        }
//
//        for (int i = 0; i < s_vv_1ring.size(); ++i) {
//            fprintf(fp, "[v <-> v] [%ld] v%d: ", s_vv_1ring[i].size(), i);
//            for (set<int>::const_iterator it = s_vv_1ring[i].begin();
//                 it != s_vv_1ring[i].end(); ++it) {
//                fprintf(fp, "%d ", *it);
//            }
//            fprintf(fp, "\n");
//        }
//
//        fclose(fp);
//    }
}

void OneRingMeshNeighborhood::computeFaceNeighbors(const Mesh *mesh) {
    m_ff_1ring.resize(mesh->getNumFaces());

    for (int f = 0; f < mesh->getNumFaces(); ++f) {
        int a = mesh->getFaces()[f].a;
        int b = mesh->getFaces()[f].b;
        int c = mesh->getFaces()[f].c;

        for (size_t j = 0; j < m_vf_1ring[a].size(); ++j) {
            m_ff_1ring[f].insert(m_vf_1ring[a][j]);
        }

        for (size_t j = 0; j < m_vf_1ring[b].size(); ++j) {
            m_ff_1ring[f].insert(m_vf_1ring[b][j]);
        }

        for (size_t j = 0; j < m_vf_1ring[c].size(); ++j) {
            m_ff_1ring[f].insert(m_vf_1ring[c][j]);
        }
    }

//    if (Config::isPrintingDebugInfo) {
//        string filename = string("DEBUG.ff_1ring.txt");
//        FILE *fp = nullptr;
//
//        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
//            fprintf(stderr, "WARNING: open DEBUG.ff_1ring.txt failed ...\n");
//            return;
//        }
//
//        for (int f = 0; f < s_ff_1ring.size(); ++f) {
//            fprintf(fp, "[f <-> f] [%ld] f%d: ", s_ff_1ring[f].size(), f);
//            for (set<int>::const_iterator it = s_ff_1ring[f].begin();
//                 it != s_ff_1ring[f].end(); ++it) {
//                fprintf(fp, "%d ", *it);
//            }
//            fprintf(fp, "\n");
//        }
//
//        fclose(fp);
//    }
}

const std::vector<std::vector<int>>& OneRingMeshNeighborhood::getVertexToFaceNeighborhood() const {
    return m_vf_1ring;
}

const std::vector<std::set<int>>& OneRingMeshNeighborhood::getVertexToVertexNeighborhood() const {
    return m_vv_1ring;
}

const std::vector<std::set<int>>& OneRingMeshNeighborhood::getFaceToFaceNeighborhood() const {
    return m_ff_1ring;
}


// TwoRingMeshNeighborhood implementations
TwoRingMeshNeighborhood::TwoRingMeshNeighborhood(const MeshNeighborhood *smallerNeighborhood) :
        OneRingMeshNeighborhood() {
    m_vf_1ring = smallerNeighborhood->getVertexToFaceNeighborhood();
    m_vv_1ring = smallerNeighborhood->getVertexToVertexNeighborhood();
    m_ff_1ring = smallerNeighborhood->getFaceToFaceNeighborhood();
    m_hasVertexToFace1RingComputed = true;
    m_hasVertexToVertex1RingComputed = true;
    m_hasFaceToFace1RingComputed = true;
}

TwoRingMeshNeighborhood::~TwoRingMeshNeighborhood() {}

void TwoRingMeshNeighborhood::computeVertexNeighbors(const Mesh *mesh) {
    if (!m_hasVertexToFace1RingComputed || !m_hasVertexToVertex1RingComputed) {
        fprintf(stderr, "Error: 1-ring vertex-face, 1-ring vertex-vertex neighborhood have not yet computed!\n");
        exit(EXIT_FAILURE);
    }

    m_vv_2ring.resize(mesh->getNumVertices());
    m_vf_2ring.resize(mesh->getNumVertices());

    // build 2-ring VERTEX <-> FACE mapping
    m_vf_2ring = m_vf_1ring;
    for (int i = 0; i < mesh->getNumFaces(); ++i) {
        // 3 vertices of current face
        int a = mesh->getFaces()[i].a;
        int b = mesh->getFaces()[i].b;
        int c = mesh->getFaces()[i].c;

        // neighboring vertices of 3 vertices
        const std::set<int> &neib1 = m_vv_1ring[a];
        const std::set<int> &neib2 = m_vv_1ring[b];
        const std::set<int> &neib3 = m_vv_1ring[c];

        // go over each neighboring vertices set, check if current face is in
        // the neighboring faces list
        std::set<int>::const_iterator it;

        for (it = neib1.begin(); it != neib1.end(); ++it) {
            std::vector<int> &neib_faces = m_vf_1ring[*it];
            if (find(neib_faces.begin(), neib_faces.end(), i) == neib_faces.end()) {
                // current face is not in the list, add it
                m_vf_2ring[*it].push_back(i);
            }
        }

        for (it = neib2.begin(); it != neib2.end(); ++it) {
            std::vector<int> &neib_faces = m_vf_1ring[*it];
            if (find(neib_faces.begin(), neib_faces.end(), i) == neib_faces.end()) {
                // current face is not in the list, add it
                m_vf_2ring[*it].push_back(i);
            }
        }

        for (it = neib3.begin(); it != neib3.end(); ++it) {
            std::vector<int> &neib_faces = m_vf_1ring[*it];
            if (find(neib_faces.begin(), neib_faces.end(), i) == neib_faces.end()) {
                // current face is not in the list, add it
                m_vf_2ring[*it].push_back(i);
            }
        }
    }

    // find 2-ring neighboring vertices
//    m_vv_2ring.assign(m_vv_1ring.begin(), m_vv_1ring.end());
    m_vv_2ring = m_vv_1ring;
    for (int i = 0; i < mesh->getNumVertices(); ++i) {
        for (int j = 0; j < m_vf_2ring[i].size(); ++j) {
            m_vv_2ring[i].insert(mesh->getFaces()[m_vf_2ring[i][j]].a);
            m_vv_2ring[i].insert(mesh->getFaces()[m_vf_2ring[i][j]].b);
            m_vv_2ring[i].insert(mesh->getFaces()[m_vf_2ring[i][j]].c);
        }
    }

//    if (Config::isPrintingDebugInfo) {
//        string filename = "DEBUG.vf_2ring.txt";
//        FILE *fp = nullptr;
//
//        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
//            fprintf(stderr, "WARNING: open DEBUG.vf_2ring.txt failed.\n");
//            return;
//        }
//
//        for (int i = 0; i < s_vf_2ring.size(); ++i) {
//            fprintf(fp, "[v <-> f] [%ld] v%d: ", s_vf_2ring[i].size(), i);
//            for (int j = 0; j < s_vf_2ring[i].size(); ++j) {
//                fprintf(fp, "%d ", s_vf_2ring[i][j]);
//            }
//            fprintf(fp, "\n");
//        }
//
//        fclose(fp);
//        fp = nullptr;
//        filename = "DEBUG.vv_2ring.txt";
//
//        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
//            fprintf(stderr, "WARNING: open DEBUG.vv_2ring.txt failed ...");
//            return;
//        }
//
//        for (int i = 0; i < s_vv_2ring.size(); ++i) {
//            fprintf(fp, "[v <-> v] [%ld] v%d: ", s_vv_2ring[i].size(), i);
//            for (set<int>::const_iterator it = s_vv_2ring[i].begin();
//                 it != s_vv_2ring[i].end(); ++it) {
//                fprintf(fp, "%d ", *it);
//            }
//            fprintf(fp, "\n");
//        }
//
//        fclose(fp);
//    }
}

void TwoRingMeshNeighborhood::computeFaceNeighbors(const Mesh *mesh) {
    if (!m_hasFaceToFace1RingComputed) {
        fprintf(stderr, "Error: 1-ring face-face has not yet computed!\n");
        exit(EXIT_FAILURE);
    }

    m_ff_2ring.resize(mesh->getNumFaces());

    for (int f = 0; f < mesh->getNumFaces(); ++f) {
        int a = mesh->getFaces()[f].a;
        int b = mesh->getFaces()[f].b;
        int c = mesh->getFaces()[f].c;

        for (int j = 0; j < m_vf_2ring[a].size(); ++j) {
            m_ff_2ring[f].insert(m_vf_2ring[a][j]);
        }

        for (int j = 0; j < m_vf_2ring[b].size(); ++j) {
            m_ff_2ring[f].insert(m_vf_2ring[b][j]);
        }

        for (int j = 0; j < m_vf_2ring[c].size(); ++j) {
            m_ff_2ring[f].insert(m_vf_2ring[c][j]);
        }
    }

//    if (Config::isPrintingDebugInfo) {
//        string filename = "DEBUG.ff_2ring.txt";
//        FILE *fp = nullptr;
//
//        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
//            fprintf(stderr, "WARNINIG: open DEBUG.ff_2ring.txt failed ...");
//            return;
//        }
//
//        for (int f = 0; f < s_ff_2ring.size(); ++f) {
//            fprintf(fp, "[f <-> f] [%ld] f%d: ", s_ff_2ring[f].size(), f);
//            for (set<int>::const_iterator it = s_ff_2ring[f].begin();
//                 it != s_ff_2ring[f].end(); ++it) {
//                fprintf(fp, "%d ", *it);
//            }
//            fprintf(fp, "\n");
//        }
//
//        fclose(fp);
//    }
}

const std::vector<std::vector<int>>& TwoRingMeshNeighborhood::getVertexToFaceNeighborhood() const {
    return m_vf_2ring;
}

const std::vector<std::set<int>>& TwoRingMeshNeighborhood::getVertexToVertexNeighborhood() const {
    return m_vv_2ring;
}

const std::vector<std::set<int>>& TwoRingMeshNeighborhood::getFaceToFaceNeighborhood() const {
    return m_ff_2ring;
}
