#include "Geometry.h"

#include <utility>


GPU_Geometry::GPU_Geometry()
	: vao()
	, vertBuffer(0, 3, GL_FLOAT)
	, colsBuffer(1, 3, GL_FLOAT)
	, normalsBuffer(2, 3, GL_FLOAT)
	, indexBuffer()
{}


void GPU_Geometry::setVerts(const std::vector<glm::vec3>& verts) {
	vertBuffer.uploadData(sizeof(glm::vec3) * verts.size(), verts.data(), GL_STATIC_DRAW);
}


void GPU_Geometry::setCols(const std::vector<glm::vec3>& cols) {
	colsBuffer.uploadData(sizeof(glm::vec3) * cols.size(), cols.data(), GL_STATIC_DRAW);
}

void GPU_Geometry::setNormals(const std::vector<glm::vec3>& norms) {
	normalsBuffer.uploadData(sizeof(glm::vec3) * norms.size(), norms.data(), GL_STATIC_DRAW);
}

void GPU_Geometry::setIndexes(const std::vector<GLuint>& indexes) {
	indexBuffer.uploadData(sizeof(GLuint) * indexes.size(), indexes.data(), GL_STATIC_DRAW);
}

