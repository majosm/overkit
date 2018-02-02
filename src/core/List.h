// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_LIST_INCLUDED
#define OVK_CORE_LIST_INCLUDED

#include "Global.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct t_list_entry {
  void *data;
  struct t_list_entry *prev;
  struct t_list_entry *next;
} t_list_entry;

typedef struct {
  t_list_entry *front;
  t_list_entry *back;
  int size;
} t_list;

static inline t_list_entry *ListBegin(t_list *List) { return List->front; }
static inline t_list_entry *ListEnd(t_list *List) { return NULL; }
static inline const t_list_entry *ListBeginC(const t_list *List) { return List->front; }
static inline const t_list_entry *ListEndC(const t_list *List) { return NULL; }

static inline t_list_entry *ListRBegin(t_list *List) { return List->back; }
static inline t_list_entry *ListREnd(t_list *List) { return NULL; }
static inline const t_list_entry *ListRBeginC(const t_list *List) { return List->back; }
static inline const t_list_entry *ListREndC(const t_list *List) { return NULL; }

static inline t_list_entry *ListFront(t_list *List) { return List->front; }
static inline t_list_entry *ListBack(t_list *List) { return List->back; }
static inline const t_list_entry *ListFrontC(t_list *List) { return List->front; }
static inline const t_list_entry *ListBackC(t_list *List) { return List->back; }

static inline int ListSize(const t_list *List) { return List->size; }

static inline t_list_entry *ListPrev(t_list_entry *Entry) { return Entry->prev; }
static inline t_list_entry *ListNext(t_list_entry *Entry) { return Entry->next; }
static inline const t_list_entry *ListPrevC(const t_list_entry *Entry) { return Entry->prev; }
static inline const t_list_entry *ListNextC(const t_list_entry *Entry) { return Entry->next; }

static inline void *ListData(t_list_entry *Entry) { return Entry->data; }
static inline const void *ListDataC(const t_list_entry *Entry) { return Entry->data; }

static inline void ListSetData(t_list_entry *Entry, void *Data) { Entry->data = Data; }

static inline void ListCreate(t_list **List_) {

  *List_ = malloc(sizeof(t_list));
  t_list *List = *List_;

  List->front = NULL;
  List->back = NULL;
  List->size = 0;

}

static inline void ListDestroy(t_list **List_) {

  t_list *List = *List_;

  OVK_DEBUG_ASSERT(List != NULL, "Invalid list pointer.");

  t_list_entry *Entry = ListBegin(List);

  while (Entry != ListEnd(List)) {
    t_list_entry *NextEntry = ListNext(Entry);
    free(Entry);
    Entry = NextEntry;
  }

  free(*List_);
  *List_ = NULL;

}

static inline t_list_entry *ListInsertBefore(t_list *List, t_list_entry *Entry, void *Data) {

  t_list_entry *NewEntry = malloc(sizeof(t_list_entry));

  NewEntry->data = Data;

  if (Entry) {
    NewEntry->prev = Entry->prev;
    NewEntry->next = Entry;
    Entry->prev = NewEntry;
  } else {
    NewEntry->prev = NULL;
    NewEntry->next = NULL;
  }

  ++List->size;

  if (!NewEntry->prev) {
    List->front = NewEntry;
  }
  if (!NewEntry->next) {
    List->back = NewEntry;
  }

  return NewEntry;

}

static inline t_list_entry *ListInsertAfter(t_list *List, t_list_entry *Entry, void *Data) {

  t_list_entry *NewEntry = malloc(sizeof(t_list_entry));

  NewEntry->data = Data;

  if (Entry) {
    NewEntry->prev = Entry;
    NewEntry->next = Entry->next;
    Entry->next = NewEntry;
  } else {
    NewEntry->prev = NULL;
    NewEntry->next = NULL;
  }

  ++List->size;

  if (!NewEntry->prev) {
    List->front = NewEntry;
  }
  if (!NewEntry->next) {
    List->back = NewEntry;
  }

  return NewEntry;

}

static inline void *ListRemove(t_list *List, t_list_entry **Entry_) {

  OVK_DEBUG_ASSERT(Entry_ != NULL, "Invalid list entry pointer.");

  t_list_entry *Entry = *Entry_;

  if (Entry->prev) {
    Entry->prev->next = Entry->next;
  }
  if (Entry->next) {
    Entry->next->prev = Entry->prev;
  }

  if (!Entry->prev) {
    List->front = Entry->next;
  }
  if (!Entry->next) {
    List->back = Entry->prev;
  }

  void *Data = Entry->data;

  free(Entry);

  *Entry_ = NULL;

  --List->size;

  return Data;

}

static inline t_list_entry *ListPushFront(t_list *List, void *Data) {

  return ListInsertBefore(List, List->front, Data);

}

static inline t_list_entry *ListPushBack(t_list *List, void *Data) {

  return ListInsertAfter(List, List->back, Data);

}

static inline void *ListPopFront(t_list *List) {

  // Don't want ListRemove call to nullify List->front
  t_list_entry *Front = List->front;

  return ListRemove(List, &Front);

}

static inline void *ListPopBack(t_list *List) {

  // Don't want ListRemove call to nullify List->back
  t_list_entry *Back = List->back;

  return ListRemove(List, &Back);

}

#ifdef __cplusplus
}
#endif

#endif
