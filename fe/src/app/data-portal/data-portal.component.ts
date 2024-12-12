import {HttpClient} from '@angular/common/http';
import { Component, AfterViewInit, ViewChild, inject } from '@angular/core';
import {MatPaginator, MatPaginatorModule} from '@angular/material/paginator';
import {MatSort, MatSortModule, SortDirection} from '@angular/material/sort';
import {MatProgressSpinnerModule} from '@angular/material/progress-spinner';
import {MatTableModule} from '@angular/material/table';
import {merge, Observable, of as observableOf} from 'rxjs';
import {catchError, map, startWith, switchMap} from 'rxjs/operators';
import {ElasticService} from "../services/elastic.service";

@Component({
  selector: 'app-data-portal',
  standalone: true,
  imports: [MatProgressSpinnerModule, MatTableModule, MatSortModule, MatPaginatorModule],
  templateUrl: './data-portal.component.html',
  styleUrl: './data-portal.component.css'
})
export class DataPortalComponent implements AfterViewInit {
  private _httpClient = inject(HttpClient);

  displayedColumns: string[] = ['urn', 'title', 'sequenceType', 'geneName', 'geneCategory', 'publicationYear',
    'numVariants'];
  elasticService: ElasticService | null | undefined;
  data: MaveDBMetadata[] = [];

  resultsLength = 0;
  isLoadingResults = true;
  isRateLimitReached = false;

  @ViewChild(MatPaginator) paginator: MatPaginator | undefined;
  @ViewChild(MatSort) sort: MatSort | undefined;

  ngAfterViewInit() {
    this.elasticService = new ElasticService(this._httpClient);

    // If the user changes the sort order, reset back to the first page.
    // @ts-ignore
    this.sort.sortChange.subscribe(() => (this.paginator.pageIndex = 0));

    // @ts-ignore
    merge(this.sort.sortChange, this.paginator.page)
      .pipe(
        startWith({}),
        switchMap(() => {
          this.isLoadingResults = true;
          return this.elasticService!.getData(
            // this.sort.active,
            // this.sort.direction,
            // this.paginator.pageIndex,
          ).pipe(catchError(() => observableOf(null)));
        }),
        map(data => {
          // Flip flag to show that loading has finished.
          this.isLoadingResults = false;
          this.isRateLimitReached = data === null;

          if (data === null) {
            return [];
          }

          // Only refresh the result length if there is new data. In case of rate
          // limit errors, we do not want to reset the paginator to zero, as that
          // would prevent users from re-triggering requests.
          this.resultsLength = data.total;
          return data.results;
        }),
      )
      .subscribe(data => (this.data = data));
  }

}

// TODO: move to another file?
export interface MaveDBMetadata {
  urn: string;
  title: string;
  sequenceType: string;
  geneName: string;
  geneCategory: string;
  publicationYear: number;
  numVariants: number;
}
