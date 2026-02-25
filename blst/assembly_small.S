#if defined(__x86_64) || defined(__x86_64__)
# if defined(__ELF__)
#  include "elf/add_mod_384-x86_64.s"
#  include "elf/add_mod_384x384-x86_64.s"
#  include "elf/add_mod_256-x86_64.s"
#  if defined(__ADX__) || defined(__BLST_PORTABLE__)
#   include "elf/ctx_inverse_mod_384-x86_64.s"
#  endif
#  if !defined(__ADX__) || defined(__BLST_PORTABLE__)
#   include "elf/ctq_inverse_mod_384-x86_64.s"
#  endif
#  include "elf/ct_inverse_mod_256-x86_64.s"
#  include "elf/div3w-x86_64.s"
#  include "elf/ct_is_square_mod_384-x86_64.s"
#  if defined(__ADX__) || defined(__BLST_PORTABLE__)
#   include "elf/mulx_mont_384-x86_64.s"
#   include "elf/mulx_mont_256-x86_64.s"
#  endif
#  if !defined(__ADX__) || defined(__BLST_PORTABLE__)
#   include "elf/mulq_mont_384-x86_64.s"
#   include "elf/mulq_mont_256-x86_64.s"
#  endif
# endif
#endif

