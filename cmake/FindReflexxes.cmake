if(NOT TARGET Reflexxes::Reflexxes)
  if(NOT ${Reflexxes_ROOT_DIR} STREQUAL "")
    if(APPLE)
      set(Reflexxes_LIB_DIR "${Reflexxes_ROOT_DIR}/MacOS/x64/release/lib/")
    elseif(UNIX)
      set(Reflexxes_LIB_DIR "${Reflexxes_ROOT_DIR}/Linux/x64/release/lib/shared")
    endif()

    set(Reflexxes_INCLUDE_DIR "${Reflexxes_ROOT_DIR}/include")
  endif()

  find_library(Reflexxes ${REFLEXXES_TYPE} PATHS ${Reflexxes_LIB_DIR})

  if(Reflexxes)
    add_library(Reflexxes::Reflexxes INTERFACE IMPORTED)
    target_include_directories(Reflexxes::Reflexxes INTERFACE ${Reflexxes_INCLUDE_DIR})
    target_link_libraries(Reflexxes::Reflexxes INTERFACE ${Reflexxes})
  endif()
endif()
