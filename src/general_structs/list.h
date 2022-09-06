// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef LIST_H
#define LIST_H

#include <stdlib.h>
#include <string.h>

/*! @brief List structure. Contains the number of items and the array of pointers to these items.
 *
 * @paragraph Details
 * The n_items is information about the number of items saves as well as allocated.
 * Memory is reallocated each time a new item is added into a list making this structure very 
 * inefficient for larger amounts of data.
 */
typedef struct list {
    size_t n_items;
    char **items;
} list_t;


/*! @brief Creates an empty list. The list must later be destroyed using list_destroy().
 * 
 * @return Pointer to list structure. If the list cannot be created, returns NULL.
 */
list_t *list_create(void);


/*! @brief Deallocates memory for the list.
 *
 * @paragraph Details
 * The memory for each individual item is also deallocated.
 * 
 * @param list      pointer to list to destroy
 */
void list_destroy(list_t *list);


/*! @brief Adds item to the end of the list.
 *
 * @paragraph Details
 * Note that list_t DOES NOT take responsibility for any dynamically allocated items.
 * If you append an item into a list, you must still deallocate the original item.
 * In other words, what is appended to the list is a copy of the original item.
 * 
 * @param list      list to change
 * @param item      string to add to the list
 * 
 * @return 0 if successful, else 1.
 */
int list_append(list_t **list, const char *item);


/*! @brief Gets an item at target index of the list. Checks boundaries.
 *
 * @paragraph Details
 * Note that you do not have to (and should not) deallocate the returned item
 * as it will be deallocated when the list is destroyed.
 * 
 * @param list      list to search
 * @param index     index of the item to return
 * 
 * @return String at the target index. If the index is out of bounds, returns NULL.
 */
char *list_get(const list_t *list, const size_t index);


/*! @brief Checks whether target string is in a list. Returns its index.
 * 
 * @paragraph Details
 * If the item is present multiple times in the list, returns the index of the first occurrence.
 *
 * @param list      list to search
 * @param item      string to search for
 * 
 * @return Non-zero if the item is present in the list. Else negative number.
 */
int list_index(const list_t *list, const char *item);

#endif /* LIST_H */