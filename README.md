# general-image-masking
General image masking approach to deal with cosmetic and structural problems in image surveys.

La idea general es tratar de mejorar las detecciones producidas en los bordes de las imagenes al utilizar
SExtractor para detectar objetos. 

Las opciones que queremos ofrecer hasta ahora:

1- Enmascarar los bordes de las imagenes.

2- Agregar bordes sinteticos a las imagenes.

3- Unir varias imagenes

Estas opciones son formas de trabajar con objetos o fuentes no astronomicas en las imagenes, de forma previa
a ser utilizadas para identificacion de objetos.

Notar que no trataremos aqui los problemas cosmeticos del CCD debidos a saturaciones o "bleddings" por estrellas saturadas.

Explicacion física del efecto que se produce sobre las estrellas saturadas en las imagenes:
https://mariotonin.me/2022/07/diffraction-patterns-in-astronomical-imaging/
