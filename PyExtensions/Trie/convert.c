#include "include/convert.h"
#include <string.h>

char convert(char c)
{ //fare stessa cosa con la funzione sotto: aggiungere case minuscoli
	switch (c)
	{
	case 'A':
		c = 0x1;
		break;
	case 'C':
		c = 0x2;
		break;
	case 'G':
		c = 0x4;
		break;
	case 'T':
		c = 0x8;
		break;
	case 'N':
		c = 0x0F;
		break;
	case '_':
		c = 0x0F;
		break;
	case 'R':
		c = 0x05;
		break;
	case 'Y':
		c = 0x0A;
		break;
	case 'S':
		c = 0x06;
		break;
	case 'W':
		c = 0x09;
		break;
	case 'K':
		c = 0x0C;
		break;
	case 'M':
		c = 0x03;
		break;
	case 'B':
		c = 0x0E;
		break;
	case 'D':
		c = 0x0D;
		break;
	case 'H':
		c = 0x0B;
		break;
	case 'V':
		c = 0x07;
		break;
	default:
		c = 0x0;
		break;
	}
	return c;
}

unsigned char convert_pam(char *c)
{
	unsigned char s;
	if (c == "AGG")
	{
		s = 0x1;
	}
	else if (c == "AGG")
	{
		s = 0x2;
	}
	else if (c == "TGG")
	{
		s = 0x3;
	}
	else if (c == "CGG")
	{
		s = 0x4;
	}
	else
	{
		s = 0x0;
	}
	free(c);
	return s;
}

char inverse(char c)
{
	switch (c)
	{
	case 0x1:
		c = 'A'; break;
	case 0x2:
		c = 'C'; break;
	case 0x4:
		c = 'G'; break;
	case 0x8:
		c = 'T'; break;
	case 0x0F:
		c = 'N'; break;
	case 0x05:
		c = 'R'; break;
	case 0x0A:
		c = 'Y'; break;
	case 0x06:
		c = 'S'; break;
	case 0x09:
		c = 'W'; break;
	case 0x0C:
		c = 'K'; break;
	case 0x03:
		c = 'M'; break;
	case 0x0E:
		c = 'B'; break;
	case 0x0D:
		c = 'D'; break;
	case 0x0B:
		c = 'H'; break;
	case 0x07:
		c = 'V'; break;
	default:
		c = 'N';
		break;
		return c;
	}
}

char* convert_dna_seq(char* seq)
{	
	int seqlen;
	if (strlen(seq) % 2 != 0) return NULL;
	char* res = (char *)malloc(strlen(seq)/2);
    for (int i =0;i<strlen(seq);i+=2){
        res[i/2] = (convert(seq[i]) << 4) + convert(seq[i+1]);
    }
	free(seq);
	return res;
}