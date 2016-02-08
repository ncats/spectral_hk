
#define __inchi_private_h__
#include "_inchi.h"

/*
 * these are from http://cccbdb.nist.gov/
 */
void
bond_energy_length (const edge_t *e,
                    double *d /* energy */, double *r /* length */)
{
  const vertex_t *u = e->u->atom->atno < e->v->atom->atno
    ? e->u : e->v;
  const vertex_t *v = _inchi_edge_other (e, u);

  /*
   * d is kJ/mol
   * r is pm
   */
  if (d != 0) *d = 0;
  if (r != 0) *r = 1000;

  switch (u->atom->atno)
    {
    case 1:
      switch (v->atom->atno)
        {
        case 1: /* H-H */
          if (d != 0) *d = 432;
          if (r != 0) *r = 74;
          break;

        case 5: /* H-B */
          if (d != 0) *d = 389;
          if (r != 0) *r = 119;
          break;

        case 6: /* H-C */
          if (d != 0) *d = 411;
          if (r != 0) *r = 109;
          break;

        case 7: /* H-N */
          if (d != 0) *d = 386;
          if (r != 0) *r = 101;
          break;

        case 8: /* H-O */
          if (d != 0) *d = 459;
          if (r != 0) *r = 96;
          break;

        case 9: /* H-F */
          if (d != 0) *d = 565;
          if (r != 0) *r = 92;
          break;

        case 14: /* H-Si */
          if (d != 0) *d = 318;
          if (r != 0) *r = 148;
          break;

        case 15: /* H-P */
          if (d != 0) *d = 322;
          if (r != 0) *r = 144;
          break;

        case 16: /* H-S */
          if (d != 0) *d = 363;
          if (r != 0) *r = 134;
          break;

        case 17: /* H-Cl */
          if (d != 0) *d = 428;
          if (r != 0) *r = 127;
          break;

        case 35: /* H-Br */
          if (d != 0) *d = 362;
          if (r != 0) *r = 141;
          break;
          
        case 53: /* H-I */
          if (d != 0) *d = 295;
          if (r != 0) *r = 161;
        }
      break;

    case 5:
      switch (v->atom->atno)
        {
        case 5: /* B-B */
          if (d != 0) *d = 293;
          if (r != 0) *r = 170.2;
          break;

        case 6: /* B-C */
          if (d != 0) *d = 536;
          if (r != 0) *r = 149.1;
          break;

        case 9:
          if (d != 0) *d = 613;
          if (r != 0) *r = 130.7;
          break;

        case 17:
          if (d != 0) *d = 456;
          if (r != 0) *r = 175;
          break;

        case 35:
          if (d != 0) *d = 377;
          if (r != 0) *r = 188.8;
          break;

        default:
          fprintf (stderr, "** Error: No energy/length entries for %s-%s! **\n",
                   u->atom->symbol, v->atom->symbol);
        }
      break;

    case 6:
      switch (v->atom->atno)
        {
        case 6:
          if (e->order == 1)
            {
              if (d != 0) *d = 346;
              if (r != 0) *r = 154;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 602;
              if (r != 0) *r = 134;
            }
          else if (e->order == 3)
            {
              if (d != 0) *d = 835;
              if (r != 0) *r = 120;
            }
          else
            fprintf (stderr, "** Error: Unknown bond order "
                     "%d between C and C **\n", e->order);
          break;

        case 7:
          if (e->order == 1)
            {
              if (d != 0) *d = 305;
              if (r != 0) *r = 147;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 615;
              if (r != 0) *r = 129;
            }
          else if (e->order == 3)
            {
              if (d != 0) *d = 887;
              if (r != 0) *r = 116;
            }
          else
            fprintf (stderr, "** Error: Unknown bond order "
                     "%d between C and N **\n", e->order);
          break;

        case 8:
          if (e->order == 1)
            {
              if (d != 0) *d = 358;
              if (r != 0) *r = 143;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 799;
              if (r != 0) *r = 120;
            }
          else if (e->order == 3)
            {
              if (d != 0) *d = 1072;
              if (r != 0) *r = 113;
            }
          else
            fprintf (stderr, "** Error: Unknown bond order "
                     "%d between C and O **\n", e->order);
          break;

        case 9: /* C-F */
          if (d != 0) *d = 485;
          if (r != 0) *r = 135;
          break;

        case 14: /* C-Si */
          if (d != 0) *d = 318;
          if (r != 0) *r = 185;
          break;

        case 15: /* C-P */
          if (d != 0) *d = 264;
          if (r != 0) *r = 184;
          break;

        case 16: /* C-S */
          if (e->order == 1)
            {
              if (d != 0) *d = 272;
              if (r != 0) *r = 182;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 573;
              if (r != 0) *r = 160;
            }
          break;
          
        case 17: /* C-Cl */
          if (d != 0) *d = 327;
          if (r != 0) *r = 177;
          break;

        case 35: /* C-Br */
          if (d != 0) *d = 285;
          if (r != 0) *r = 194;

        case 53: /* C-I */
          if (d != 0) *d = 213;
          if (r != 0) *r = 214;
        }
      break;

    case 7:
      switch (v->atom->atno)
        {
        case 7: /* N*N */
          if (e->order == 1)
            {
              if (d != 0) *d = 167;
              if (r != 0) *r = 145;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 418;
              if (r != 0) *r = 125;
            }
          else if (e->order == 3)
            {
              if (d != 0) *d = 942;
              if (r != 0) *r = 110;
            }
          break;

        case 8: /* N*O */
          if (e->order == 1)
            {
              if (d != 0) *d = 201;
              if (r != 0) *r = 140;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 607;
              if (r != 0) *r = 121;
            }
          break;

        case 9: /* N-F */
          if (d != 0) *d = 283;
          if (r != 0) *r = 136;
          break;

        case 14: /* N-Si */
          if (d != 0) *d = 355;
          if (r != 0) *r = 157.19;
          break;

        case 16: /* N=S */
          if (r != 0) *r = 149.7;
          break;
          
        case 17: /* N-Cl */
          if (d != 0) *d = 313;
          if (r != 0) *r = 175;
          break;
        }
      break;

    case 8:
      switch (v->atom->atno)
        {
        case 8:
          if (e->order == 1)
            {
              if (d != 0) *d = 142;
              if (r != 0) *r = 148;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 494;
              if (r != 0) *r = 121;
            }
          break;

        case 9:
          if (d != 0) *d = 190;
          if (r != 0) *r = 142;
          break;

        case 14:
          if (d != 0) *d = 452;
          if (r != 0) *r = 163;
          break;

        case 15:
          if (e->order == 1)
            {
              if (d != 0) *d = 335;
              if (r != 0) *r = 163;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 544;
              if (r != 0) *r = 150;
            }
          break;

        case 16:
          if (e->order == 1)
            {
              if (r != 0) *r = 157.4;
            }
          else if (e->order == 2)
            {
              if (d != 0) *d = 522;
              if (r != 0) *r = 143;
            }
          break;

        case 53:
          if (d != 0) *d = 201;
          break;
        }
      break;

    case 9:
      switch (v->atom->atno)
        {
        case 9:
          if (d != 0) *d = 155;
          if (r != 0) *r = 142;
          break;
          
        case 14:
          if (d != 0) *d = 565;
          if (r != 0) *r = 160;

        case 15:
          if (d != 0) *d = 490;
          if (r != 0) *d = 154;
          break;

        case 16:
          if (d != 0) *d = 284;
          if (r != 0) *r = 156;
        }
      break;

    case 14:
      switch (v->atom->atno)
        {
        case 14:
          if (d != 0) *d = 222;
          if (r != 0) *r = 233;
          break;
          
        case 16:
          if (d != 0) *d = 293;
          if (r != 0) *r = 200;
          break;

        case 17:
          if (d != 0) *d = 381;
          if (r != 0) *r = 202;
          break;

        case 35:
          if (d != 0) *d = 310;
          if (r != 0) *r = 215;
          break;

        case 53:
          if (d != 0) *d = 234;
          if (r != 0) *r = 243;
          break;
        }
      break;
      
    case 15:
      switch (v->atom->atno)
        {
        case 15:
          if (d != 0) *d = 201;
          if (r != 0) *r = 221;
          break;
          
        case 16:
          if (e->order == 2)
            {
              if (d != 0) *d = 335;
              if (r != 0) *r = 186;
            }
          break;

        case 17:
          if (d != 0) *d = 326;
          if (r != 0) *r = 203;
          break;
          
        case 35:
          if (d != 0) *d = 264;
          break;

        case 53:
          if (d != 0) *d = 184;
        }
      break;

    case 16:
      switch (v->atom->atno)
        {
        case 16:
          if (e->order == 2)
            {
              if (d != 0) *d = 425;
              if (r != 0) *r = 149;
            }
          break;

        case 17:
          if (d != 0) *d = 255;
          if (r != 0) *r = 207;
          break;
        }
      break;

    case 17:
      if (v->atom->atno == 17)
        {
          if (d != 0) *d = 240;
          if (r != 0) *r = 199;
        }
      else if (v->atom->atno == 53)
        {
          if (d != 0) *d = 208;
          if (r != 0) *r = 232;
        }
      break;

    case 35:
      if (v->atom->atno == 35)
        {
          if (d != 0) *d = 190;
          if (r != 0) *r = 228;
        }
      else if (v->atom->atno == 53)
        {
          if (d != 0) *d = 175;
        }
      break;

    case 53:
      if (v->atom->atno == 53)
        {
          if (d != 0) *d = 148;
          if (r != 0) *r = 267;
        }
      break;
    }
}
