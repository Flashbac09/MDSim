#ifndef _PRODUCT_H_
#define _PRODUCT_H_

// extern void all_verlet_list_standard();
// extern void all_verlet_list_extended();
//憋用啦

extern double get_distance(int i, int j);

/*
extern void all_verlet_list_standard_sorted();
extern void all_verlet_list_extended_sorted();
*/
//这两个不要了，在old分支里可以用
/*
extern void specific_verlet_list_standard(int n);
extern void specific_verlet_list_extended(int n);
extern void specific_verlet_list_standard_sorted(int n);
extern void specific_verlet_list_extended_sorted(int n);
*/
// specific的不想重写了,没用
extern void nlist();

#endif