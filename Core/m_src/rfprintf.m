function rfprintf(st)

  global st__;
  
  erase= repmat('%c', [1 numel(st__)]);
  fprintf(erase, 8*ones(1, numel(st__)));
  fprintf(st);
  st__ = st;
