PUBLISHED_PDF = ../../../tutorials/tutorial_2/slides.pdf

slides.pdf: slides.tex
	xelatex -shell-escape slides.tex

slides.tex: slides.org
	emacs --script ./org-files-to-tex.el slides.org

force:
	xelatex -8bit -shell-escape slides.tex

manual:
	xelatex -shell-escape slides.tex

publish: slides.pdf
	mkdir -p $(dir $(PUBLISHED_PDF))
	cp slides.pdf $(PUBLISHED_PDF)

watch:
	@# Re-run automatically if slides.tex or slides.org changes.
	latexmk --xelatex --shell-escape -pvc -view=none slides.tex

clean:
	rm -f *.tex slides.pdf *.tex~ *.nav *.log *.out *.snm *.toc *.aux *.vrb *.fls *.fdb_latexmk *.synctex.gz
	rm -rf _minted-*

final:
	make clean
	make slides.pdf
	make force
	make publish

help:
	@echo "make slides.pdf -- compile org to tex to pdf"
	@echo "make slides.tex -- update slides.tex from slides.org"
	@echo "make publish    -- copy slides.pdf to the root tutorials folder"
	@echo "make watch      -- update slides.pdf on changes"
	@echo "make clean      -- delete derived files"
	@echo "make final      -- rebuild from a clean working directory"
	@echo "make help       -- show this command overview"
