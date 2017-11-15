// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ORDERED_MAP_INCLUDED
#define OVK_CORE_ORDERED_MAP_INCLUDED

#include "Debug.h"
#include "Global.h"

typedef struct t_ordered_map_entry {
  int key;
  void *data;
  struct t_ordered_map_entry *next;
} t_ordered_map_entry;

typedef struct {
  t_ordered_map_entry *head;
  int size;
} t_ordered_map;

static inline t_ordered_map_entry *OMBegin(t_ordered_map *Map) { return Map->head; }
static inline t_ordered_map_entry *OMEnd(t_ordered_map *Map) { return NULL; }
static inline const t_ordered_map_entry *OMBeginC(const t_ordered_map *Map) { return Map->head; }
static inline const t_ordered_map_entry *OMEndC(const t_ordered_map *Map) { return NULL; }

static inline int OMSize(const t_ordered_map *Map) { return Map->size; }

static inline t_ordered_map_entry *OMNext(t_ordered_map_entry *Entry) { return Entry->next; }
static inline const t_ordered_map_entry *OMNextC(const t_ordered_map_entry *Entry) { return Entry->next; }

static inline int OMKey(const t_ordered_map_entry *Entry) { return Entry->key; }

static inline void *OMData(t_ordered_map_entry *Entry) { return Entry->data; }
static inline const void *OMDataC(const t_ordered_map_entry *Entry) { return Entry->data; }

static inline void OMCreate(t_ordered_map **Map_) {

  *Map_ = malloc(sizeof(t_ordered_map));
  t_ordered_map *Map = *Map_;

  Map->head = NULL;
  Map->size = 0;

}

static inline void OMDestroy(t_ordered_map **Map_) {

  t_ordered_map *Map = *Map_;

  OVK_DEBUG_ASSERT(Map != NULL, "Invalid ordered map pointer.");

  t_ordered_map_entry *Entry = OMBegin(Map);

  while (Entry != OMEnd(Map)) {
    t_ordered_map_entry *NextEntry = OMNext(Entry);
    free(Entry);
    Entry = NextEntry;
  }

  free(*Map_);
  *Map_ = NULL;

}

static inline int OMNextAvailableKey(const t_ordered_map *Map) {

  const t_ordered_map_entry *Entry = OMBeginC(Map);

  if (Entry == OMEndC(Map) || Entry->key > 0) {
    return 0;
  }

  const t_ordered_map_entry *PrevEntry = Entry;
  Entry = OMNextC(PrevEntry);
  while (Entry != OMEndC(Map)) {
    if (Entry->key - PrevEntry->key > 1) {
      return PrevEntry->key + 1;
    }
    PrevEntry = Entry;
    Entry = OMNextC(PrevEntry);
  }

  return PrevEntry->key + 1;

}

static inline t_ordered_map_entry *OMFind(t_ordered_map *Map, int Key) {

  t_ordered_map_entry *Entry = OMBegin(Map);

  while (Entry != OMEnd(Map) && Entry->key < Key) {
    Entry = OMNext(Entry);
  }

  if (Entry != OMEnd(Map) && Entry->key == Key) {
    return Entry;
  } else {
    return OMEnd(Map);
  }

}

static inline const t_ordered_map_entry *OMFindC(const t_ordered_map *Map, int Key) {

  const t_ordered_map_entry *Entry = OMBeginC(Map);

  while (Entry != OMEndC(Map) && Entry->key < Key) {
    Entry = OMNextC(Entry);
  }

  if (Entry != OMEndC(Map) && Entry->key == Key) {
    return Entry;
  } else {
    return OMEndC(Map);
  }

}

static inline bool OMExists(const t_ordered_map *Map, int Key) {

  const t_ordered_map_entry *Entry = OMFindC(Map, Key);

  return Entry != OMEndC(Map);

}

static inline t_ordered_map_entry *OMInsert(t_ordered_map *Map, int Key, void *Data) {

  t_ordered_map_entry *NewEntry = malloc(sizeof(t_ordered_map_entry));
  NewEntry->key = Key;
  NewEntry->data = Data;
  NewEntry->next = NULL;

  t_ordered_map_entry *Entry = OMBegin(Map);

  if (Entry == OMEnd(Map) || Entry->key > Key) {
    Map->head = NewEntry;
    NewEntry->next = Entry;
  } else {
    t_ordered_map_entry *PrevEntry = Entry;
    Entry = OMNext(PrevEntry);
    while (Entry != OMEnd(Map) && Entry->key < Key) {
      PrevEntry = Entry;
      Entry = OMNext(PrevEntry);
    }
    OVK_DEBUG_ASSERT(Entry == OMEnd(Map) || Entry->key != Key, "Entry with specified key already "
      "exists.");
    NewEntry->next = PrevEntry->next;
    PrevEntry->next = NewEntry;
  }

  ++Map->size;

  return NewEntry;

}

static inline void *OMRemove(t_ordered_map *Map, int Key) {

  t_ordered_map_entry *Entry = OMBegin(Map);

  OVK_DEBUG_ASSERT(Entry != OMEnd(Map) && Entry->key <= Key, "Entry with specified key does not "
    "exist.");

  if (Entry->key == Key) {
    Map->head = Entry->next;
  } else {
    t_ordered_map_entry *PrevEntry = Entry;
    Entry = OMNext(PrevEntry);
    while (Entry != OMEnd(Map) && Entry->key < Key) {
      PrevEntry = Entry;
      Entry = OMNext(PrevEntry);
    }
    OVK_DEBUG_ASSERT(Entry->key == Key, "Entry with specified key does not exist.");
    PrevEntry->next = Entry->next;
  }

  void *Data = Entry->data;

  free(Entry);

  --Map->size;

  return Data;

}

#endif
