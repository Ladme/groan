// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

/* Implementation of a dictionary (hash table). */

#ifndef DICT_H
#define DICT_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>


/*! @brief Entry to a dictionary (key->value pair). */
typedef struct entry {
    char *key;
    void *value;
} dict_entry_t;


/*! @brief Dictionary (hash table) */
typedef struct dict {
    size_t allocated;       // the number of entries for which memory has been allocated
    size_t available;       // the number of positions for entries that are available
    dict_entry_t *entries;
} dict_t;


/*! @brief Creates an empty dictionary and allocates memory for it.
 * 
 * @return Pointer to dictionary structure. If the dictionary cannot be created, returns NULL.
 */
dict_t *dict_create(void);


/*! @brief Frees memory allocated to dictionary.
 *
 * @paragraph Notes
 * If dict is NULL, this function does nothing.
 *
 * @param dict          pointer to dictionary to be destroyed
 */
void dict_destroy(dict_t *dict);


/*! @brief Add target key->value pair into the dictionary
 *
 * @paragraph Multiple identical keys
 * If multiple identical keys are added into the dictionary, only the last
 * value associated with this key will be preseverd.
 *
 * @param dict          pointer to dictionary to work with
 * @param key           key to set
 * @param value         value to set (pointer!)
 * @param valsize       size of the value in bytes
 * 
 * @return Zero if the key->value pair was successfully added. Else non-zero.
 */
int dict_set(dict_t *dict, const char *key, const void *value, const unsigned valsize);


/*! @brief Gets the value of target key
 *
 * @param dict          pointer to dictionary to work with
 * @param key           key to search for
 * 
 * @return Void pointer to the value. Null if the entry has not been found.
 */ 
void *dict_get(const dict_t *dict, const char *key);

/*! @brief Gets all keys of target dictionary
 *
 * @paragraph Note on deallocation
 * The resulting 'keys' array must be later deallocated using simple free().
 * (The individual strings will be deallocated when the associated dictionary is destroyed.)
 * 
 * @param dict          pointer to dictionary to work with
 * @param keys          pointer to an array of strings
 * 
 * @return Number of keys in the dictionary
 */ 
size_t dict_keys(dict_t *dict, char ***keys);


#endif /* DICT_H */