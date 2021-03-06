#include "util/wt/incl/util.h"
#include "multal.h"

struct	slot	{
		int	id;
		struct slot	*more, *less;
		}
	;
typedef struct	slot	Slot;

struct	best	{
		int	score, m, l;
		struct best	*more, *less;
		}
	;
typedef struct	best	Best;

typedef struct	{
		int	id;
		int	count;
		int	*slot;
		}
	Pept;
Pept	*pept_list;
int	n_pept = 0, n_prot = 0;

struct node	{
		int	id;
		int	count;
		struct	node	*more, *less;
		Slot	*slot;
		}
	;
typedef struct	node	Node;
Node	*key;

struct pair	{
		int	id;
		char	count;
		struct	pair	*more, *less;
		}
	;
typedef struct	pair	Pair;
Pair	*gate;
int	ngates,
	free_pairs = 0;

typedef struct	{
		int	id;
		unsigned char	count;
		}
	Prot;
Prot	*pair_list;
int	n_pair = 0;

int	base[PEPLEN],
	n_alloc, p_alloc = 1000,
	jump, joff;

int	soft_level = 0, diag_up, min_pept, max_hold;

int	*nth;

int	get_nth();
Slot	*get_slot();
Best	*get_best();
Node	*get_node();
Pair	*get_pair();
char	*decode();
char	soften();

init_slots (nprots,pep_len,soft)
int	nprots, pep_len, soft;
{	int	peplen = abs(pep_len),
		i, j, power = 1;
	printf("\nCalculating peptide table based on length %d\n", peplen);
	if (soft) {
		printf("peptide softened to level %d\n", soft);
		soft_level = soft;
		if (pep_len<0) {
			printf("(default softening ignored)\n");
			pep_len = -pep_len;
		}
	}
	if (pep_len<0) {
		soft_level = peplen-4;
		printf("peptide softened by default to level %d\n", soft_level);
	}
	for (i=peplen; i>0; i--) {
		base[i-1] = power;
		power = power*NACIDS;
	}
	nth = (int*)malloc(sizeof(int)*(nprots+1));
	TEST(nth)
	key = (Node*)malloc(sizeof(Node)*(nprots+1));
	TEST(key)
	for (i=1; i<=nprots; i++) {
	{	Node	*k = key+i;
		nth[i] = min_pept;
		k->id = i;
		k->slot = 0;
		k->count = 0;
		k->more = k->less = 0;
		}
	}
}

char	soften(res,soft)
char	res;
int	soft;
{
	res = UPPER(res);
	if (res=='B') return 'D';
	if (res=='Z') return 'E';
	if (soft==0) return res;
	if (res=='S') return 'T';
	if (res=='L') return 'I';
	if (soft==1) return res;
	if (res=='D') return 'E';
	if (res=='N') return 'Q';
	if (soft==2) return res;
	if (res=='V') return 'I';
	if (res=='R') return 'K';
	if (soft==3) return res;
	if (res=='Y') return 'F';
	if (res=='A') return 'G';
	return res;
}

find_slot (open,pep_len,pep_in,prot)
char	*pep_in;
int	open, pep_len, prot;
{	int	peplen = abs(pep_len), 
		coff = 'A'-1,
		hole, id = 0, i;
	Node	*p = key+prot;
	char	peptide[PEPLEN];
	for (i=0; i<peplen; i++) {
		peptide[i] = soften(pep_in[i],soft_level);
		id += (peptide[i]-coff)*base[i];
	}
	peptide[peplen] = '\0';
	if (p->slot) get_slot(p,p->slot,id);
	else	p->slot = get_slot(p,0,id);
}

Slot *get_slot (node, slot, new_id)
Node	*node;
Slot	*slot;
int	new_id;
{
	Slot 	*n = slot;
	if (slot) 
	{	int	slot_id = slot->id;
		if(new_id == slot_id) {
			return n;
		}
		if (new_id < slot_id) {
			if (slot->less) get_slot(node,slot->less,new_id);
			else	slot->less = get_slot(node,0,new_id);
		} else {
			if (slot->more) get_slot(node,slot->more,new_id);
			else	slot->more = get_slot(node,0,new_id);
		}
	} else {
		n = (Slot*)malloc(sizeof(Slot));
		TEST(n)
		n->id = new_id;
		n->less = n->more = 0;
		(node->count)++;
	}
	return n;
}

int	get_nth (slot,n)
Best	*slot;
int	*n;
{ int	s = min_pept;
	if (!slot) return;
	if (slot->more) s = get_nth(slot->more,n);
	(*n)++;
	if (*n==max_hold) return slot->score;
	if (slot->less) s = get_nth(slot->less,n);
	return s;
}

free_best (slot)
Best	*slot;
{
	if (!slot) return;
	if (slot->less) free_best(slot->less);
	if (slot->more) free_best(slot->more);
	free (slot);
}

in_slot (id, jd)
int	id;
{
	int	i, j, n;
	Pept	*pid = pept_list+id,
		*pjd = pept_list+jd;
	i = j = n = 0;
	while (i < pid->count && j < pjd->count) {
		if (pid->slot[i]  < pjd->slot[j]) i++;
		if (pid->slot[i]  > pjd->slot[j]) j++;
		if (pid->slot[i] == pjd->slot[j]) {
			i++; j++; n++;
		}
	}
	return n;
}

Best *get_best (slot, n, more, kept)
Best	*slot;
int	n, more, *kept;
{
	Best	*s = slot;
	if (slot) {
		if (n <= slot->score) {
			more += slot->m + 1;
			if (slot->less) {
				slot->less = get_best(slot->less,n,more,kept);
				slot->l = (slot->less)->l + (slot->less)->m + 1;
			} else {
				slot->less = get_best(0,n,more,kept);
				if (slot->less) slot->l = 1;
			}
		} else {
			if (slot->more) {
				slot->more = get_best(slot->more,n,more,kept);
				slot->m = (slot->more)->m + (slot->more)->l + 1;
			} else {
				slot->more = get_best(0,n,more,kept);
				if (slot->more) slot->m = 1;
			}
			if (slot->more && more+slot->m >= max_hold) {
				s = slot->more;
				if (slot->less) free_best(slot->less);
				free (slot);
			}
		}
	} else {
		if (more >= max_hold) return 0;
		s = (Best*)malloc(sizeof(Best));
		TEST(s)
		s->score = n;
		s->less = s->more = 0;
		s->l = s->m = 0;
		*kept = 1;
	}
	return s;
}

count (nprots,aspan,bspan,place,minp,hold,diag)
int	nprots, aspan, bspan, minp, hold, diag;
int	*place;
{
	int	n, i, j, k, left, deep, topi, topj, n_pairs,
		npairs = (nprots*nprots-nprots)/2,
		pept_max, pair_limit = 16000000;
	Best	**best;
	best = (Best**)malloc(sizeof(Best*)*(nprots+1));
	for (i=1; i<=nprots; i++) best[i] = 0;
	min_pept = minp;
	max_hold = hold;
	if (diag < 0) diag_up = min_pept;
		else  diag_up = diag;
	ngates = nprots;
	jump = npairs/ngates;
	left = npairs-jump*ngates;
	j = joff = (jump+left+1)/2;
	if (!place[0] || free_pairs)
		gate = (Pair*)malloc(sizeof(Pair)*(ngates+1));
	for (i=0; i<ngates; i++) {
	{	Pair	*g = gate+i;
			g->id = j;
			g->count = 0;
			g->more = g->less = 0;
			j += jump;
		}
	}
	printf("\n\n ****** CYCLE %d ******\n", place[0]);
	printf("Minimum score held = %d\n", min_pept);
	printf("Top %d scores held per sequence\n", max_hold);
	printf("Sequences in the range %d<Dij<%d scored\n",aspan,bspan);
	printf("Adjacent sequences scored with bonus of %d\n", diag_up);
	n_pairs = 0;
	for (i=1; i<nprots; i++)
	{ int	pid = place[i],
		pjd = place[i+1],
		id = pack(pid,pjd,nprots),
		way = min(max(0,(id-joff+1)/jump),nprots-1),
		n = in_slot(pid,pjd);
		get_pair(gate+way,id,-n,place[0]);
		npairs++;
	}
	for (i=1; i<=(nprots-aspan); i++)
	{ int	pid = place[i];
		for (j=i+aspan; j<=min(nprots,i+bspan); j++) 
		{ int	pjd = place[j],
			pdif = abs(i-j),
			id = pack(pid,pjd,nprots),
			way = min(max(0,(id-joff+1)/jump),nprots-1),
			n = in_slot(pid,pjd);
			if (n_pair == pair_limit) break;
			if (n<min_pept) continue;
			topi = topj = 0;
			if (n>nth[pid]) {
				if (!best[pid])
					best[pid] = get_best(0,n,0,&topi);
				else best[pid] = get_best(best[pid],n,0,&topi);
				deep = 0;
				if (topi) nth[pid] = get_nth(best[pid],&deep);
			}
			if (n>nth[pjd]) {
				if (!best[pjd])
					best[pjd] = get_best(0,n,0,&topj);
				else best[pjd] = get_best(best[pjd],n,0,&topj);
				deep = 0;
				if (topj) nth[pjd] = get_nth(best[pjd],&deep);
			}
			if (topi || topj) {
				get_pair(gate+way,id,n,place[0]);
				n_pairs++;
			}
			if (n_pair == pair_limit) break;
 		}
		free_best(best[pid]);
		if (n_pair == pair_limit) break;
	}
	printf("%d pairs collected\n",n_pairs);
	if (n_pair == pair_limit) printf("%d pair slots filled\n", pair_limit);
	return pair_out();
}

Pair *get_pair (pair, id, score, again)
Pair	*pair;
int	id, score;
int	again;
{
	if (pair)
	{	int	pair_id = pair->id;
		if (id == pair_id) {
			if (pair->count && again) return;
			if (pair->count+score < BYTE)
				pair->count += score;
			else	pair->count  = BYTE;
			return;
		}	 
		if (id < pair_id) {
			if (pair->less) get_pair(pair->less,id,score,again);
			else	pair->less = get_pair(0,id,score,0);
		} else {
			if (pair->more) get_pair(pair->more,id,score,again); 	
			else	pair->more = get_pair(0,id,score,0);
		}
	} 
	else 
	{	Pair	*n;
		n = (Pair*)malloc(sizeof(Pair));
		TEST(n)
		n->id = id;
		n->less = n->more = 0;
		n->count = score;
		n_pair++;
		return n;
	}
}

pscore (pid, pjd, nprots)
int	pid, pjd;
{	int	id = pack(pid,pjd,nprots),
		way = (id-joff+1)/jump,
		count;
	way = min(max(0,way),nprots-1);
	count = find_pair(gate+way,id);
	return count;
 }

find_pair (pair, id)
Pair	*pair;
int	id;
{	int	count, pair_id = pair->id;
	if (id == pair_id) return (int)pair->count;
	count = 0;
	if (id < pair_id) {
		if (pair->less) count = find_pair(pair->less,id);
	} else {
		if (pair->more) count = find_pair(pair->more,id); 	
	}
	return count;
}

pair_out()
{
int	i, j;
	n_pair = 0;
	n_alloc = p_alloc;
	pair_list = (Prot*)malloc(sizeof(Prot)*n_alloc);
	TEST(pair_list)
	for (i=0; i<ngates; i++) print_pair(gate+i);
	if (free_pairs) free(gate);
	PRINTi(n_pair) NL
	return n_pair;
}

print_pair (pair)
Pair	*pair;
{	
	if (pair->less) print_pair(pair->less);
	if (abs(pair->count) >= min_pept)
	{	Prot	*p;
		if (n_pair==n_alloc) {
			n_alloc += p_alloc;
			pair_list = (Prot*)
				realloc(pair_list,n_alloc*sizeof(Prot));
			TEST(pair_list)
		}
		pair_list[n_pair].id = pair->id;
		if (pair->count < 0) {
			pair_list[n_pair].count = diag_up - pair->count;
			pair->count = -(pair->count);
		} else  pair_list[n_pair].count = pair->count;
		n_pair++;
	}
	if (pair->more) print_pair(pair->more);
	if (free_pairs) free(pair);
	return;
}

pept_out(nprots)
int	nprots;
{
int	i, j;
	pept_list = (Pept*)malloc(sizeof(Pept)*(nprots+1));
	TEST(pept_list)
	for (i=1; i<=nprots; i++)
	{	Node	*node = key+i;
		Pept	*pep = pept_list+i;
		pep->id = node->id;
		pep->count = node->count;
		pep->slot = (int*)malloc(pep->count*sizeof(int));
		TEST(pep->slot)
		n_prot = 0;
		if (node->slot) print_slot(node->slot,pep);
		else printf("*NB* no slot for prot = %d\n",i);
	}
	free(key);
}

print_slot (slot,pep)
Slot	*slot;
Pept	*pep;
{
	if (slot->less) print_slot(slot->less,pep);
	if (slot->id) pep->slot[n_prot++] = slot->id;
	if (slot->more) print_slot(slot->more,pep);
	free(slot);
	return;
}

char	*decode (n,peplen)
int	n, peplen;
{
	char 	pep[PEPLEN];
	int	i, c, rem = n;
	for (i=0; i<peplen; i++) {
		c = rem/base[i];
		pep[i]='A'+c-1;
		rem -= c*base[i];
	}
	pep[peplen] = '\0';
	return pep;
}	
reorder (pairs)
int	*pairs;
{
	int	i;
	int	*score;
	int	*place;
printf("reorder ");
	score = (int*)  malloc(n_pair*sizeof(int));
		TEST(score)
	place = (int*)malloc(n_pair*sizeof(int));
		TEST(place)
printf(" malloc OK");
	for (i=0; i<n_pair; i++) score[i] = pair_list[i].count;
	sort(0, score, 0, place, n_pair, 1);
	free(score);
	for (i=0; i<n_pair; i++) pairs[i] = pair_list[place[i]].id;
	free(pair_list);
	free(place);
printf(" done\n");
	return;
}
