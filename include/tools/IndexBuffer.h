#pragma once

#include "GLHandles.h"

#include <glad/glad.h>


class IndexBuffer {

public:
	IndexBuffer();

	// Because we're using the IndexBufferHandle to do RAII for the buffer for us
	// and our other types are trivial or provide their own RAII
	// we don't have to provide any specialized functions here. Rule of zero
	//
	// https://en.cppreference.com/w/cpp/language/rule_of_three
	// https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#Rc-zero

	// Public interface
	void bind() const { glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferID); }
	void uploadData(GLsizeiptr size, const void* data, GLenum usage);

private:
	IndexBufferHandle bufferID;
};