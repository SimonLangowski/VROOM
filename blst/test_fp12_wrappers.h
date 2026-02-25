/*
 * Wrapper header to make static functions from fp12_tower.c accessible for testing
 * Include this after including fp12_tower.c to access internal functions
 */

#ifndef TEST_FP12_WRAPPERS_H
#define TEST_FP12_WRAPPERS_H

// We'll need to modify fp12_tower.c to export these, or create wrapper functions
// For now, we'll create wrapper functions that call the internal logic

// Wrapper for sqr_fp4 - we'll need to implement this or make the static function public
// Since sqr_fp4 is static, we'll create a test version that duplicates the logic
// or we can test it indirectly through functions that use it

#endif

