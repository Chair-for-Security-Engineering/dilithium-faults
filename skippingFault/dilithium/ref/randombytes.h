#ifndef RANDOMBYTES_H
#define RANDOMBYTES_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#include <unistd.h>
#include <sys/syscall.h>
#endif

#include <stddef.h>
#include <stdint.h>

void randombytes(uint8_t *out, size_t outlen);

#endif
