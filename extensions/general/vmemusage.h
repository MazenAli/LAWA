#ifndef EXTENSIONS_GENERAL_VMEMUSAGE_H
#define EXTENSIONS_GENERAL_VMEMUSAGE_H 1

namespace extensions
{

int
parseLine(char* line);

/* Note: this value is in KB! */
int
vmemusage(void);

double
kbtogb(int kb);

} // namespace extensions

#endif // EXTENSIONS_GENERAL_VMEMUSAGE_H
