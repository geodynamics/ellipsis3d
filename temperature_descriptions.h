/*

Copyright (C) 1995 The GeoFramework Consortium

This file is part of Ellipsis3D.

Ellipsis3D is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License, version 2,
as published by the Free Software Foundation.

Ellipsis3D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Author:
  Louis Moresi <louis.moresi@sci.monash.edu>

*/



struct TEMP_OPT {
    void (* update_temperature)();

  struct Rect Trects;
  struct Circ Tcircs;
  struct Harm Tharms;
  struct RectBc Trectbcs;
  struct CircBc Tcircbcs;
  struct HarmBc Tharmbcs;
  struct PolyBc Tpolybcs;

  struct Rect Tintz_on;
  struct Rect Tintz_off;
    
    
} temperature;

