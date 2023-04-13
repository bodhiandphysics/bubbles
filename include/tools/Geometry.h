#pragma once

//------------------------------------------------------------------------------
// This file contains simple classes for storing geomtery on the CPU and the GPU
// Later assignments will require you to expand these classes or create your own
// similar classes with the needed functionality
//------------------------------------------------------------------------------

#include "VertexArray.h"
#include "VertexBuffer.h"
#include "IndexBuffer.h"

#include <glad/glad.h>
#include <glm/glm.hpp>

#include <vector>


// VAO and two VBOs for storing vertices and texture coordinates, respectively
class GPU_Geometry {

public:
	GPU_Geometry();

	// Public interface
	void bind() { vao.bind(); }
	void readydraw() {vao.bind(); indexBuffer.bind();}

	void setVerts(const std::vector<glm::vec3>& verts);
	void setIndexes(const std::vector<GLuint>& indexes);
	void setCols(const std::vector<glm::vec3>& cols);
	void setNormals(const std::vector<glm::vec3>& norms);

private:
	// note: due to how OpenGL works, vao needs to be
	// defined and initialized before the vertex buffers
	VertexArray vao;

	VertexBuffer vertBuffer;
	VertexBuffer normalsBuffer;
	VertexBuffer colsBuffer;
	IndexBuffer  indexBuffer;
};
