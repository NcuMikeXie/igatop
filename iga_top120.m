%% A 120 LINE ISOGEOMETRIC TOPOLOGY OPTIMIZATION CODE BASED ON BEZIER EXTRACTION %%
function iga_top120(nelx,nely,volfrac,penal,rmin)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-3;
nu = 0.3;
%% PREPARE IGA
noPtsX = nelx + 2;
noPtsY = nely + 2;
noCtrPts = noPtsX * noPtsY;
noElems = nelx * nely;
noDofs = 2*noCtrPts;
noPtsX_B = 2*nelx + 1;
noPtsY_B = 2*nely + 1;
noCtrPts_B = noPtsX_B * noPtsY_B;
A11 = [144,-24,-40,24,-36,-28,-8,-20,-12;-24,128,-24,-36,32,-36,-20,0,-20;-40,-24,144,-28,-36,24,-12,-20,-8;
        24,-36,-28,112,-8,-24,24,-36,-28;-36,32,-36,-8,96,-8,-36,32,-36;-28,-36,24,-24,-8,112,-28,-36,24;
        -8,-20,-12,24,-36,-28,144,-24,-40;-20,0,-20,-36,32,-36,-24,128,-24;-12,-20,-8,-28,-36,24,-40,-24,144];
A12 = [45,-30,-15,30,-20,-10,15,-10,-5;30,0,-30,20,0,-20,10,0,-10;15,30,-45,10,20,-30,5,10,-15;
      -30,20,10,0,0,0,30,-20,-10;-20,0,20,0,0,0,20,0,-20;-10,-20,30,0,0,0,10,20,-30;
      -15,10,5,-30,20,10,-45,30,15;-10,0,10,-20,0,20,-30,0,30;-5,-10,15,-10,-20,30,-15,-30,45];
B11 = [-48,-24,-8,24,12,4,24,12,4;-24,-32,-24,12,16,12,12,16,12;-8,-24,-48,4,12,24,4,12,24;
        24,12,4,-48,-24,-8,24,12,4;12,16,12,-24,-32,-24,12,16,12;4,12,24,-8,-24,-48,4,12,24;
        24,12,4,24,12,4,-48,-24,-8;12,16,12,12,16,12,-24,-32,-24;4,12,24,4,12,24,-8,-24,-48];
B12 = [45,90,45,-90,-20,-10,-45,-10,-5;-90,0,90,20,0,-20,10,0,-10;-45,-90,-45,10,20,90,5,10,45;
       90,20,10,0,0,0,-90,-20,-10;-20,0,20,0,0,0,20,0,-20;-10,-20,-90,0,0,0,10,20,90;
       45,10,5,90,20,10,-45,-90,-45;-10,0,10,-20,0,20,90,0,-90;-5,-10,-45,-10,-20,-90,45,90,45];
C14 = full(sparse([1 4 7 2 5 8 3 6 9],1:9,ones(1,9)));
KE = E0/(1-nu^2)/360*([A11 A12;A12' C14'*A11*C14]+nu*[B11 B12; B12' C14'*B11*C14]);
Cxi = repmat([0.5,0,0;0.5,1,0.5;0,0,0.5],1,1,nelx);
Cxi(:,1,1) = [1;0;0];
Cxi(:,3,nelx) = [0;0;1];
Cet = repmat([0.5,0,0;0.5,1,0.5;0,0,0.5],1,1,nely);
Cet(:,1,1) = [1;0;0];
Cet(:,3,nely) = [0;0;1];
CxiG = sparse(noPtsX,noPtsX_B);
CxiG(noPtsX*([1 2:2:noPtsX_B-1 noPtsX_B]-1)+(1:noPtsX))=1;
CxiG(noPtsX*([3:2:noPtsX_B-2 3:2:noPtsX_B-2]-1)+[2:noPtsX-2,3:noPtsX-1])=0.5;
CetG = sparse(noPtsY,noPtsY_B);
CetG(noPtsY*([1 2:2:noPtsY_B-1 noPtsY_B]-1)+(1:noPtsY))=1;
CetG(noPtsY*([3:2:noPtsY_B-2 3:2:noPtsY_B-2]-1)+[2:noPtsY-2,3:noPtsY-1])=0.5;
elConnU_B = bsxfun(@plus,(1:2:nelx*2-1)',0:2);
elConnV_B = bsxfun(@plus,(1:2:nely*2-1)',0:2);
element_B = (reshape(permute(repmat(elConnV_B,1,1,nelx,3),[3,1,4,2]),noElems,9)-1)*noPtsX_B+repmat(elConnU_B,nely,3);
edofMat_B = [element_B noCtrPts_B+element_B];
elConnU = bsxfun(@plus,(1:nelx)',0:2);
elConnV = bsxfun(@plus,(1:nely)',0:2);
element = (reshape(permute(repmat(elConnV,1,1,nelx,3),[3,1,4,2]),noElems,9)-1)*noPtsX+repmat(elConnU,nely,3);
edofMat = [element noCtrPts+element];
iK = reshape(kron(edofMat,ones(18,1))',324*noElems,1);
jK = reshape(kron(edofMat,ones(1,18))',324*noElems,1);
clear elConnU_B elConnV_B element_B elConnU elConnV element edofMat
%% BOUNDARY OF HALF MBB BEAM
F = sparse(noCtrPts+noPtsX*(noPtsY-1)+1,1,-1,noDofs,1);
U = zeros(noDofs,1);
fixeddofs = union(1:noPtsX:noCtrPts, noCtrPts+noPtsX);
alldofs = 1:noDofs;
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(noElems*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nely
  for j1 = 1:nelx
    e1 = (i1-1)*nelx+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nely)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nelx)
        e2 = (i2-1)*nelx+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,nelx,nely);
loop = 0;
change = 1;
while change > 0.01
  loop = loop + 1;
  %% IGA
  sK = zeros(324*noElems,1,1);
  k = 0;
  for j = 1:nely
      for i = 1:nelx
          k = k+1;
          CK = [kron(Cet(:,:,j),Cxi(:,:,i)),zeros(9,9);zeros(9,9),kron(Cet(:,:,j),Cxi(:,:,i))];
          sK((k-1)*324+1:k*324) = reshape(CK*KE*CK',324,1)*(Emin+x(i,j).^penal*(E0-Emin));
      end
  end
  K = sparse(iK,jK,sK);
  K=(K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  U_B = full([kron(CetG,CxiG)'*U(1:noCtrPts);kron(CetG,CxiG)'*U(noCtrPts+(1:noCtrPts))]);
  %% OBJECTIVE FUNCTION CALCULATION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U_B(edofMat_B)*KE).*U_B(edofMat_B),2),nelx,nely);
  c = sum(sum((Emin+x.^penal*(E0-Emin)).*ce));
  dc = -penal*(E0-Emin)*x.^(penal-1).*ce;
  dv = ones(nelx,nely);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    if sum(xnew(:)) > volfrac*noElems, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(x(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-flipud(x')); caxis([0 1]); axis equal; axis off; drawnow;
end
% =========================================================================
% === This code was written by X XIE and A Yang, School of Advanced     ===
% === Manufacturing, Nanchang University, Nanchang, CHINA               ===
% === ----------------------------------------------------------------- ===
% === Please send your suggestions and comments to: xiexd020@ncu.edu.cn ===
% === ----------------------------------------------------------------- ===