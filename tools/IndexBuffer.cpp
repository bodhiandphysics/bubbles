#include "IndexBuffer.h"

#include <utility>


IndexBuffer::IndexBuffer()
	: bufferID{}
{}


void IndexBuffer::uploadData(GLsizeiptr size, const void* data, GLenum usage) {
	if (size > 0) {
		bind();
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, data, usage);
	}

}