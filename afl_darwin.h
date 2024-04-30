#include <string.h>
#include <stdbool.h>
#include <stdint.h>

#ifndef DARWIN_H
#define DARWIN_H

typedef unsigned int uint;

void DARWIN_init(uint32_t, unsigned);

int DARWIN_SelectOperator(uint32_t);

uint32_t DARWIN_NotifyFeedback(uint32_t, unsigned);

uint32_t DARWIN_get_parent_repr(uint32_t);

#endif
