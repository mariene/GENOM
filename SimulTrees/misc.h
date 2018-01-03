#ifndef _MISC_H_
#define _MISC_H_



struct elt {
	int id;
	struct elt *prev;
	struct elt *next;
	int pos;
};

#endif
