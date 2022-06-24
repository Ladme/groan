// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

/* Implementation of a dictionary (hash table). */

#include "dict.h"

/* CONSTANTS */

/*! \brief Hashing constant for hash_key function */
static const unsigned long FNV_OFFSET = 14695981039346656037UL;
/*! \brief Hashing constant for hash_key function */
static const unsigned long FNV_PRIME = 1099511628211UL;

/*! \brief The number of key->item pairs for which memory is initially allocated. */
static const size_t DICT_BLOCK = 64;

/******************************************************/

/* STATIC FUNCTIONS */


/*! \brief Hashing function for string. */
static uint64_t hash_key(const char *key)
{
    uint64_t hash = FNV_OFFSET;
    for (size_t i = 0; key[i]; ++i) {
        hash ^= (uint64_t)(unsigned char)(key[i]);
        hash *= FNV_PRIME;
    }

    return hash;
}


/*! \brief Allocates memory for dict entry and assigns key and value. Returns 0 if successful, else returns 1. */
static int dict_entry_create(dict_entry_t *entry, const char *key, const void *value, const unsigned valsize)
{
    entry->key = calloc(1, strlen(key) + 1);
    if (entry->key == NULL) {
        return 1;
    }
    memcpy(entry->key, key, strlen(key) + 1);

    entry->value = calloc(1, valsize);
    if (entry->value == NULL) {
        free(entry->key);
        return 1;
    }
    memcpy(entry->value, value, valsize);

    return 0;
}


/*! \brief Deallocates memory for a dictionary entry. */
static void dict_entry_destroy(dict_entry_t *entry)
{
    free(entry->key);
    free(entry->value);
}


/*! \brief Loops through an array of entries searching for suitable position.
 * 
 * \par Suitable position
 * Suitable position is either an empty position or position containing the same key as 'key'.
 * 
 * \par No suitable position
 * If there is no suitable position in the array, this function loops indefinitely.
 * 
 * \return Index of a suitable position.
 */
static size_t dict_find_suitable(
        const dict_entry_t *entries, 
        const size_t len, 
        const size_t starting_index,
        const char *key)
{
    size_t index = starting_index;
    while (entries[index].key != NULL && strcmp(entries[index].key, key) != 0) {
        ++index;
        // wrap to the start of array if the index iterator reaches the end
        if (index >= len) index = 0;
    }

    return index;
}

/******************************************************/

dict_t *dict_create(void) 
{
    // allocate memory for the dict itself
    dict_t *dict = calloc(1, sizeof(dict_t));
    if (dict == NULL) {
        return NULL;
    }

    // allocate memory for the array of entries
    dict->entries = calloc(DICT_BLOCK, sizeof(dict_entry_t));
    if (dict->entries == NULL) {
        free(dict);
        return NULL;
    }

    // set the initial size of the dictionary
    dict->allocated = DICT_BLOCK;
    dict->available = DICT_BLOCK / 2;

    return dict; 
}

void dict_destroy(dict_t *dict)
{
    if (dict == NULL) return;

    // loop through all entries and deallocate memory
    for (size_t i = 0; i < dict->allocated; ++i) {
        dict_entry_destroy(&dict->entries[i]); 
    }
    
    free(dict->entries);
    free(dict);
}

int dict_set(dict_t *dict, const char *key, const void *value, const unsigned valsize)
{
    // if the number of available positions is 0 (or below), we must allocate more memory
    // the number of available positions does not have to (and should not) correspond to
    // the number of _allocated_ positions
    if (dict->available <= 0) {
        size_t new_alloc = dict->allocated * 2;

        // allocate new array for the entries
        dict_entry_t *new_entries = calloc(new_alloc, sizeof(dict_entry_t));
        if (new_entries == NULL) {
            dict_destroy(dict);
            return 1;
        }

        // loop through the original entries and recalculate their indices
        // then copy these entries to their new positions
        for (size_t i = 0; i < dict->allocated; ++i) {
            if (dict->entries[i].key != NULL) {
                size_t new_index = hash_key(dict->entries[i].key) % (new_alloc);
                // new index may not be empty, we have to check that and change the index if needed
                new_index = dict_find_suitable(new_entries, new_alloc, new_index, dict->entries[i].key);
                memcpy(&(new_entries[new_index]), &(dict->entries[i]), sizeof(dict_entry_t));
            }
        }

        // we have to free the original entries
        free(dict->entries);
        dict->entries = new_entries;
        dict->available += (new_alloc - dict->allocated) / 2;
        dict->allocated = new_alloc;
    }

    // get usable index for the entry
    size_t index = hash_key(key) % dict->allocated;
    index = dict_find_suitable(dict->entries, dict->allocated, index, key);
    // remove any data that may be on the position
    dict_entry_destroy(&(dict->entries[index]));

    // once we find suitable position, populate the entry
    if (dict_entry_create(&(dict->entries[index]), key, value, valsize) != 0) {
        dict_destroy(dict);
        return 1;
    }

    // decrease the number of available positions
    --dict->available;

    return 0;
}

void *dict_get(dict_t *dict, const char *key)
{
    size_t index = hash_key(key) % dict->allocated;

    dict_entry_t *entry = &(dict->entries[index]);
    // loop forwards in the array of entries until we find the corresponding entry or empty position
    // once we reach an empty position, we can be sure that the entry does not exist in the dictionary
    while (entry->key != NULL && strcmp(entry->key, key)) {
        ++index;
        // wrap to the start of the array
        if (index >= dict->allocated) index = 0;
        entry = &(dict->entries[index]);
    }

    // return the value of the found entry
    // (this of course actually returns void pointer to the value)
    // this returns NULL, if no corresponding entry has been found
    return entry->value;
}