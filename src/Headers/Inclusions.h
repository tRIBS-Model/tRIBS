/*******************************************************************************
 * TIN-based Real-time Integrated Basin Simulator (tRIBS)
 * Distributed Hydrologic Model
 *
 * Copyright (c) 2025. tRIBS Developers
 *
 * See LICENSE file in the project root for full license information.
 ******************************************************************************/

/***************************************************************************
**
**  Inclusions.h: General tRIBS Include Statements
**
***************************************************************************/

#ifndef INCLUSIONS_H
#define INCLUSIONS_H


// INCLUDED LIBRARY HEADER FILES

#ifdef ALPHA_64
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <assert.h>
#include <memory.h>
#include <unistd.h>

#elif defined LINUX_32
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <memory>
#include <vector>
#include <unistd.h>

#elif defined MAC
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <memory>
#include <vector>
#include <unistd.h>

#elif defined WIN
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <assert.h>
#include <memory.h>

#else
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <assert.h>
#include <memory.h>

#endif

// INCLUDED FILES 

#include "src/Headers/Definitions.h"
#include "src/Headers/Classes.h"
#include "src/Headers/globalFns.h"
#include "src/tHydro/tHydroModel.h"
#include "src/tFlowNet/tFlowNet.h"
#include "src/tMesh/tMesh.h"
#include "src/tMeshList/tMeshList.h"
#include "src/tMeshElements/meshElements.h"
#include "src/tPtrList/tPtrList.h"
#include "src/tList/tList.h"
#include "src/tListInputData/tListInputData.h"
#include "src/tCNode/tCNode.h"
#include "src/Mathutil/mathutil.h"
#include "src/tArray/tArray.h"
#include "src/tArray/tMatrix.h"
#include "src/tInOut/tInputFile.h"
#include "src/tInOut/tOutput.h"
#include "src/tInOut/tOstream.h"
#include "src/tRasTin/tResample.h"
#include "src/tStorm/tStorm.h"
#include "src/tRasTin/tRainfall.h"
#include "src/tRasTin/tRainGauge.h"
#include "src/tRasTin/tVariant.h"
#include "src/tRasTin/tInvariant.h"
#include "src/tSimulator/tRunTimer.h"
#include "src/tSimulator/tControl.h"
#include "src/tHydro/tIntercept.h"
#include "src/tHydro/tEvapoTrans.h"
#include "src/tHydro/tHydroMet.h"
#include "src/tHydro/tHydroMetStoch.h"
#include "src/Mathutil/predicates.h"
#include "src/Mathutil/geometry.h"
#include "src/tHydro/tWaterBalance.h"
#include "src/tSimulator/tPreProcess.h"

#endif

//=========================================================================
//
//
//                   End of Inclusions.h  
//
//
//=========================================================================
