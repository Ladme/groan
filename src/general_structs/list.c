// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

/* Implementation of a simple dynamic array for storing strings.*/

#include "list.h"
//#include <stdio.h>

list_t *list_create(void)
{
    // allocate memory for the list
    list_t *list = calloc(1, sizeof(list_t));
    if (list == NULL) return NULL;

    return list;
}

void list_destroy(list_t *list)
{
    // loop throgh items of the list and deallocate them
    for (size_t i = 0; i < list->n_items; ++i) {
        free(list->items[i]);
    }

    // deallocate the entire list
    free(list->items);
    free(list);
}

int list_append(list_t **list, const char *item)
{
    // reallocate memory for the list
    char **new_items = realloc((*list)->items, ((*list)->n_items + 1) * sizeof(char *));
    if (new_items == NULL) return 1;

    (*list)->items = new_items;
    
    // allocate memory for the item pointer
    (*list)->items[(*list)->n_items] = calloc(strlen(item) + 1, 1);

    // copy item to the list
    strcpy( (*list)->items[(*list)->n_items], item);
    (*list)->n_items++;

    return 0;
}

char *list_get(const list_t *list, const size_t index)
{
    if (index >= list->n_items) return NULL;
    
    return list->items[index];
}

int list_index(const list_t *list, const char *item)
{
    for (size_t i = 0; i < list->n_items; ++i) {
        if (!strcmp(list->items[i], item)) return (int) i;
    }

    return -1;
}

/*
int main(void)
{
    list_t *list = list_create();

    char **items = (char *[]){"1", "test", "test2", "something", "blabla", "x"};

    for (size_t i = 0; i < 6; ++i) {
        printf("%s\n", items[i]);
        list_append(&list, items[i]);
    }

    for (size_t i = 0; i < list->n_items; ++i) {
        printf("Item %ld: %s\n", i, list_get(list, i));
    }

    if (list_get(list, 6) == NULL) printf("NONEXISTENT\n");
    if (list_get(list, 1000) == NULL) printf("NONEXISTENT\n");
    if (list_get(list, -1) == NULL) printf("NONEXISTENT\n");
    if (list_get(list, -187) == NULL) printf("NONEXISTENT\n");

    printf("%d ", (list_index(list, "test")));
    printf("%d ", (list_index(list, "blabla")));
    printf("%d\n", (list_index(list, "does not exist")));

    list_destroy(list);

    return 0;
}*/