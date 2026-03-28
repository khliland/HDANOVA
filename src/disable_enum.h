#if defined(__APPLE__)
# if defined(__clang__)
#  undef HAVE_ENUM_BASE_TYPE
#  define HAVE_ENUM_BASE_TYPE 0
# endif
#endif
